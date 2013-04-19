// Copyright (c) 2013 Genome Research Limited.
//
// This file is part of Sam Untangle.
//
// Sam Untangle is free software: you can redistribute it and/or modify it under the
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
    char* merged_input_name;
    char* unaccounted_header_name;
    char* unaccounted_name;
    size_t input_count;
    char** input_name;
    char** input_header_name;
    char** input_output_name;
};

typedef struct parsed_opts parsed_opts_t;

struct state {
    samFile* merged_input_file;
    bam_hdr_t* merged_input_header;
    samFile* unaccounted_file;
    bam_hdr_t* unaccounted_header;
    size_t input_count;
    FILE** input_file;
    samFile** input_output_file;
    bam_hdr_t** input_output_header;
};

typedef struct state state_t;

void cleanup_state(state_t* status);
void cleanup_opts(parsed_opts_t* opts);


parsed_opts_t* parse_args(int argc, char** argv) {
    if (argc < 3) {
        dprintf(STDERR_FILENO, "Arguments should be: fix_merge <merged.bam> <unaccounted_header.sam> <unaccounted.bam> <input1readnames.txt:input1header.sam:input1output.bam> <input2readnames.txt:input2header.sam:input2output.bam> [<inputXreadnames.txt:inputXheader.sam:inputXoutput.bam> ...]\r\n");
        return NULL;
    }

    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) return NULL;

    retval->merged_input_name = strdup(argv[1]);
    retval->unaccounted_header_name = strdup(argv[2]);
    retval->unaccounted_name = strdup(argv[3]);

    retval->input_count = argc-4;
    retval->input_name = (char**)calloc(retval->input_count,sizeof(char*));
    retval->input_header_name = (char**)calloc(retval->input_count,sizeof(char*));
    retval->input_output_name = (char**)calloc(retval->input_count,sizeof(char*));
    for (size_t i = 0; i < retval->input_count; i++) {
        char* temp = strdup(argv[i+4]);
        char* sep = temp;
        retval->input_name[i] = strsep(&sep, ":");
        retval->input_header_name[i] = strsep(&sep, ":");
        retval->input_output_name[i] = sep;
    }

    return retval;
}

state_t* init(parsed_opts_t* opts) {
    state_t* retval = malloc(sizeof(state_t));
    if (!retval) {
        dprintf(STDERR_FILENO, "Out of memory");
        return NULL;
    }

    retval->merged_input_file = sam_open(opts->merged_input_name, "r", 0);
    if (!retval->merged_input_file) {
        dprintf(STDERR_FILENO, "Could not open header file (%s)", opts->merged_input_name);
        return NULL;
    }
    retval->merged_input_header = sam_hdr_read(retval->merged_input_file);
    samFile* hdr_load = sam_open(opts->unaccounted_header_name, "r", 0);
    if (!hdr_load) {
        dprintf(STDERR_FILENO, "Could not open unaccounted header file (%s)", opts->unaccounted_header_name);
        return NULL;
    }
    retval->unaccounted_header = sam_hdr_read(hdr_load);
    sam_close(hdr_load);
    
    retval->unaccounted_file = sam_open(opts->unaccounted_name, "wb", 0);
    if (retval->unaccounted_file == NULL) {
        dprintf(STDERR_FILENO, "Could not open unaccounted output file: %s\r\n", opts->unaccounted_name);
        cleanup_state(retval);
        return NULL;
    }

    // Open files from repeated option
    retval->input_count = opts->input_count;
    retval->input_file = (FILE**)calloc(opts->input_count, sizeof(FILE*));
    retval->input_output_file = (samFile**)calloc(opts->input_count, sizeof(samFile*));
    if (!retval->input_file || !retval->input_output_file) {
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
        
        retval->input_output_file[i] = sam_open(opts->input_output_name[i], "wb", 0);
        samFile* headerFile = sam_open(opts->input_header_name[i], "r", 0);
        if (!hdr_load) {
            dprintf(STDERR_FILENO, "Could not open input header file (%s)", opts->input_header_name[i]);
            return NULL;
        }
        retval->input_output_header[i] = sam_hdr_read(headerFile);
        sam_close(headerFile);
    }

    return retval;
}

bool fix_merge(state_t* opts) {
    if (sam_hdr_write(opts->unaccounted_file, opts->unaccounted_header) != 0) {
        dprintf(STDERR_FILENO, "Could not write output file header");
        return false;
    }
    for (size_t i = 0; i < opts->input_count; i++) {
        if (sam_hdr_write(opts->input_output_file[i], opts->input_output_header[i]) != 0) {
            dprintf(STDERR_FILENO, "Could not write output file header");
            return false;
        }
    }
    
    bam1_t* file_read = bam_init1();
    // Read the first record
    if (sam_read1(opts->merged_input_file, opts->merged_input_header, file_read) < 0) {
        // Nothing more to read?  Ignore this file
        bam_destroy1(file_read);
        file_read = NULL;
    }
    
    // Fill initial read names to compare
    char **input_read_name = (char**)calloc(opts->input_count, sizeof(char*));
    size_t *input_buff_len = (size_t*)calloc(opts->input_count, sizeof(size_t));
    for (size_t i = 0; i < opts->input_count; i++) {
        if( getline(&input_read_name[i],&input_buff_len[i], opts->input_file[i]) < 0 ) {
            free(input_read_name[i]);
            input_read_name[i] = NULL;
        }
    }

    while (file_read != NULL) {
        bool found = false;
        size_t i;
        for (i = 0; i < opts->input_count; i++) {
            if ((input_read_name[i] != NULL) && !strcmp(bam_get_qname(file_read), input_read_name[i])) {
                found = true;
                // replace read name we just ate
                if( getline(&input_read_name[i],&input_buff_len[i], opts->input_file[i]) < 0 ) {
                    free(input_read_name[i]);
                    input_read_name[i] = NULL;
                }

                break;
            }
        }
        
        // Write the read out and replace it with the next one to process
        if (found) {
            // if found write to the appropriate untangled bam
            sam_write1(opts->input_output_file[i], opts->input_output_header[i], file_read);
        } else {
            // otherwise write to the unaccounted bam
            sam_write1(opts->unaccounted_file, opts->unaccounted_header, file_read);
        }
        
        if (sam_read1(opts->merged_input_file, opts->merged_input_header, file_read) < 0) {
            // Nothing more to read?  Ignore this file in future
            bam_destroy1(file_read);
            file_read = NULL;
        }
    }

    return true;
}

void cleanup_state(state_t* status) {
    sam_close(status->unaccounted_file);
    sam_close(status->merged_input_file);
    for (size_t i = 0; i < status->input_count; i++) {
        fclose(status->input_file[i]);
        sam_close(status->input_output_file[i]);
    }
    free(status->input_file);
}

void cleanup_opts(parsed_opts_t* opts) {
    free(opts->unaccounted_name);
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
