// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of fix_merge.
//
// fix_merge is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#define _GNU_SOURCE

#include <htslib/sam.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

struct parsed_opts {
    char* output_header_name;
    char* merged_input_name;
    size_t input_count;
    char** input_name;
    char** input_header_name;
    char** input_trans_name;
    char* output_name;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
    samFile* merged_input_file;
    bam_hdr_t* merged_input_header;
    size_t input_count;
    FILE** input_file;
    int32_t** input_trans;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct state state_t;

void cleanup_state(state_t* status);
void cleanup_opts(parsed_opts_t* opts);


parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        dprintf(STDERR_FILENO, "Arguments should be: fix_merge <newheader.sam> <merged.bam> <input1readnames.txt:input1header.sam:trans_tbl.txt> <input2readnames.txt:input2header.sam:trans_tbl.txt> [<inputXreadnames.txt:inputXheader.sam:trans_tbl.txt> ...] <output.bam>\r\n");
        return NULL;
    }

    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;

    retval->output_header_name = strdup(argv[1]);
    retval->merged_input_name = strdup(argv[2]);

    retval->input_count = argc-4;
    retval->input_name = (char**)calloc(retval->input_count,sizeof(char*));
    retval->input_header_name = (char**)calloc(retval->input_count,sizeof(char*));
    retval->input_trans_name = (char**)calloc(retval->input_count,sizeof(char*));
    size_t i = 0;
    for (; i < retval->input_count; i++) {
        char* temp = strdup(argv[i+3]);
        char* sep = temp;
        retval->input_name[i] = strsep(&sep, ":");
        retval->input_header_name[i] = strsep(&sep, ":");
        retval->input_trans_name[i] = sep;
    }

    retval->output_name = strdup(argv[i+3]);

    return retval;
}

int* build_translation_file(const char* trans_name, bam_hdr_t* file_header, bam_hdr_t* replace_header) {
    FILE* trans_file = fopen(trans_name, "r");
    
    int file_entries = file_header->n_targets;
    int replace_entries = replace_header->n_targets;
    
    int *trans = malloc(sizeof(int)*file_entries);
    
    char* linepointer = NULL;
    size_t read = 0;
    
    int counter = file_entries;
    
    while (!feof(trans_file) && !ferror(trans_file) && counter > 0) {
        getline(&linepointer, &read, trans_file);
        
        char* sep = linepointer;
        strsep(&sep, "\t");
        if (sep == NULL) break;
        char* two = sep;
        strsep(&two, "\t\n");

        // lookup tid of original and replacement
        int i = 0;
        for ( ; i < file_entries; i++ ) {
            char* item = file_header->target_name[i];
            if (!strcmp(item,linepointer)) { break; }
        }
        int j = 0;
        for ( ; j < replace_entries; j++ ) {
            char* item = replace_header->target_name[j];
            if (!strcmp(item,sep)) { break; }
        }
        
        trans[i] = j;
        counter--;
    }
    free(linepointer);
    
    fclose(trans_file);
    return trans;
}

state_t* init(parsed_opts_t* opts) {
    state_t* retval = malloc(sizeof(state_t));
    if (!retval) {
        dprintf(STDERR_FILENO, "Out of memory");
        return NULL;
    }

    // Now the input files
    samFile* hdr_load = sam_open(opts->output_header_name, "r", 0);
    if (!hdr_load) {
        dprintf(STDERR_FILENO, "Could not open header file (%s)", opts->output_header_name);
        return NULL;
    }
    retval->output_header = sam_hdr_read(hdr_load);
    sam_close(hdr_load);

    retval->merged_input_file = sam_open(opts->merged_input_name, "r", 0);
    if (!retval->merged_input_file) {
        dprintf(STDERR_FILENO, "Could not open header file (%s)", opts->merged_input_name);
        return NULL;
    }
    retval->merged_input_header = sam_hdr_read(retval->merged_input_file);

    retval->input_count = opts->input_count;
    
    // Open files
    retval->input_trans = (int**)calloc(opts->input_count, sizeof(int*));
    retval->input_file = (FILE**)calloc(opts->input_count, sizeof(samFile*));
    if (!retval->input_file || !retval->input_trans) {
        dprintf(STDERR_FILENO, "Out of memory");
        free(retval);
        return NULL;
    }
    for (size_t i = 0; i < opts->input_count; i++) {
        retval->input_file[i] = fopen(opts->input_name[i], "rb");
        if (retval->input_file[i] == NULL) {
            dprintf(STDERR_FILENO, "Could not open input file: %s\r\n", opts->input_name[i]);
            return NULL;
        }
        samFile* headerFile = sam_open(opts->input_header_name[i], "r", 0);
        if (!hdr_load) {
            dprintf(STDERR_FILENO, "Could not open input header file (%s)", opts->input_header_name[i]);
            return NULL;
        }

        bam_hdr_t* header = sam_hdr_read(headerFile);
        retval->input_trans[i] = build_translation_file(opts->input_trans_name[i], header, retval->output_header);
        free(header);
        sam_close(headerFile);
    }

    // Open output files
    retval->output_file = sam_open(opts->output_name, "wb", 0);
    if (retval->output_file == NULL) {
        dprintf(STDERR_FILENO, "Could not open output file: %s\r\n", opts->output_name);
        cleanup_state(retval);
        return NULL;
    }
    
    return retval;
}

bool fix_merge(state_t* opts) {
    if (sam_hdr_write(opts->output_file, opts->output_header) != 0) {
        dprintf(STDERR_FILENO, "Could not write output file header");
        return false;
    }
    
    bam1_t* file_read = bam_init1();
    // Read the first record
    if (sam_read1(opts->merged_input_file, opts->merged_input_header, file_read) < 0) {
        // Nothing more to read?  Ignore this file
        bam_destroy1(file_read);
        file_read = NULL;
    }
    
    char **input_read_name = (char**)calloc(opts->input_count, sizeof(char*));
    size_t *input_len = (size_t*)calloc(opts->input_count, sizeof(size_t));
    for (size_t i = 0; i < opts->input_count; i++) {
        if( getline(&input_read_name[i],&input_len[i], opts->input_file[i]) < 0 ) {
            free(input_read_name[i]);
            input_read_name[i] = NULL;
        }
    }

    while (file_read != NULL) {
        bool found = false;
        size_t i;
        for (i = 0; i < opts->input_count; i++) {
            if (!strcmp(bam_get_qname(file_read), input_read_name[i])) {
                found = true;
                break;
            }
        }
        
        if (found) {
            // if found transform the input file TID and MTID
            if (file_read->core.tid != -1) {
                file_read->core.tid = opts->input_trans[i][file_read->core.tid];
            }
            if (file_read->core.mtid != -1) {
                file_read->core.mtid = opts->input_trans[i][file_read->core.mtid];
            }
            if( getline(&input_read_name[i],&input_len[i], opts->input_file[i]) < 0 ) {
                free(input_read_name[i]);
                input_read_name[i] = NULL;
            }
        }
        
        // Write the read out and replace it with the next one to process
        sam_write1(opts->output_file, opts->output_header, file_read);
        if (sam_read1(opts->merged_input_file, opts->merged_input_header, file_read) < 0) {
            // Nothing more to read?  Ignore this file in future
            bam_destroy1(file_read);
            file_read = NULL;
        }
    }

    return true;
}

void cleanup_state(state_t* status) {
    sam_close(status->output_file);
    sam_close(status->merged_input_file);
    for (size_t i = 0; i < status->input_count; i++) {
        fclose(status->input_file[i]);
    }
    free(status->input_file);
}

void cleanup_opts(parsed_opts_t* opts) {
    free(opts->output_name);
    for (size_t i = 0; i < opts->input_count; i++) {
        free(opts->input_name[i]);
    }
    free(opts->input_name);
}

int main(int argc, char** argv) {

    parsed_opts_t* opts = parse_args(argc, argv);
    if (!opts ) return -1;
    state_t* status = init(opts);
    if (!status) return -1;
    
    if (!fix_merge(status)) return -1;
    
    cleanup_state(status);
    cleanup_opts(opts);
      
    return 0;
}
