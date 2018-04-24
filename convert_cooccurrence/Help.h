/*
 * Copyright 2018 Frangou Lab
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#ifndef Help_h
#define Help_h

#include <string>
#include <iostream>

inline void PrintHelp(FILE *destination)
{
    printf("\
NAME:\n\
    coocc2inter v0.1\n\
\n\
OVERVIEW:\n\
    Converts gene co-occurrence matrix into a gene interaction table.\n\
    Interactions will be ordered by decreasing absolute value of each one's\n\
    strength.\n\
\n\
USAGE:\n\
    coocc2inter <input table> [-o <output path>]\n\
\n\
OPTIONS:\n\
    -h                     - Show this message.\n\
    -v                     - Verbose output.\n\
    -f                     - Override the output file even if it already exists.\n\
\n\
EXAMPLES:\n\
    coocc2inter ~/input.csv                      - Convert co-occurrence matrix '~/input.csv' into a gene\n\
                                                   interaction table '~/input-interactions.csvc'.\n\
    coocc2inter ~/input.csv -o ~/output.csvc -v  - Convert co-occurrence matrix '~/input.csv' into a gene\n\
                                                   interaction table '~/output.csvc'. Verbose mode will\n\
                                                   enable reporting genes that don't have any interactions\n\
                                                   with other genes in the table.\n\
");
}

class ArgumentsParser {
 public:
    std::string input_file_path;
    std::string output_file_path;
    
    bool verbose_output{false};
    bool override_output{false};
    
    ArgumentsParser(int argc, const char * argv[])
    {
        if (argc == 1) {
            PrintHelp(stdout);
            return;
        }
        
        for (int i = 1; i < argc; ++i) {
            std::string argument = argv[i];
            if (argument == "-v" || argument == "--verbose") {
                verbose_output = true;
            } else if (argument == "-o" || argument == "--output") {
                i++;
                if (i == argc || argv[i][0] == '-') {
                    fprintf(stderr, "No valid output path has been entered.\n");
                    std::abort();
                }
                output_file_path.assign(argv[i]);
            } else if (argument == "-h" || argument == "--help") {
                PrintHelp(stdout);
                std::exit(0);
            } else if (argument == "-f" || argument == "--force") {
                override_output = true;
            } else if (argument[0] != '-') {
                input_file_path.assign(argument);
            } else {
                PrintHelp(stderr);
                fprintf(stderr, "Unknown flag '%s'\n", argument.c_str());
                std::exit(1);
            }
        }
    }
};

#endif /* Help_h */
