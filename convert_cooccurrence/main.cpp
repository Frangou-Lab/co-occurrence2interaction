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

#include "Help.h"
#include "CoOccurrenceReporter.hpp"

#include "../libgene/source/file/TxtFile.hpp"
#include "../libgene/source/file/sequence/SequenceFile.hpp"
#include "../libgene/source/utils/Tokenizer.hpp"
#include "../libgene/source/utils/CppUtils.hpp"
#include "../libgene/source/utils/StringUtils.hpp"
#include "../libgene/source/def/Flags.hpp"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, const char *argv[])
{
    ArgumentsParser arguments(argc, argv);
    
    if (arguments.input_file_path.empty()) {
        fprintf(stderr, "No input files provided. Terminating\n");
        return 1;
    }

    std::unique_ptr<SequenceFile> input_file;
    try {
        auto flags = std::make_unique<CommandLineFlags>();
        flags->SetSetting(Flags::kInputFormat, "txt");
        input_file = SequenceFile::FileWithName(arguments.input_file_path,
                                                flags,
                                                OpenMode::Read);
    } catch (std::runtime_error err) {
        if (!input_file) {
            fprintf(stderr,
                    "File \'%s\' couldn't be opened. Either it doesn't exist, or you don't have permissions to read it.\n",
                    arguments.input_file_path.c_str());
            return 1;
        }
    }
    
    size_t dot_position = arguments.input_file_path.rfind('.');
    std::string extension = utils::GetExtension(arguments.input_file_path);
    
    if (extension != "csvc" && extension != "tsvc") {
        // Force the output file to be interpreted as a 'columns-defined' file.
        extension += 'c';
    }
    
    if (arguments.output_file_path.empty())
        arguments.output_file_path = arguments.input_file_path.substr(0, dot_position) + "-interactions" + "." + extension;
    
    char delimiter = ',';
    if (extension == "tsv" || extension == "tsvc") {
        delimiter = '\t';
    }
    
    std::unique_ptr<SequenceFile> out_file;
    if (!arguments.override_output) {
        FILE *test_out_file = fopen(arguments.output_file_path.c_str(), "wx");
        if (!test_out_file) {
            fprintf(stdout,
                    "File \'%s\' already exists. Do you wish to override it? [Y/n] ",
                    arguments.output_file_path.c_str());
            char response;
            std::cin >> response;
            if (std::tolower(response) != 'y') {
                fprintf(stderr, "Skipping file \'%s\'\n", arguments.input_file_path.c_str());
                return 1;
            }
        } else {
            fclose(test_out_file);
        }
    }
    
    auto flags = std::make_unique<CommandLineFlags>();
    flags->SetSetting(Flags::kOutputFormat, "txt");
    if (!(out_file = SequenceFile::FileWithName(arguments.output_file_path,
                                                flags,
                                                OpenMode::Write))) {
        fprintf(stderr, "Couldn't open the output file \'%s\'\n",
                arguments.output_file_path.c_str());
        return 1;
    }
    
    // Write header to output file.
    SequenceRecord header;
    header.seq = "Gene Symbol,Interaction Gene Symbol";
    out_file->Write(header);
    
    SequenceRecord input_record = input_file->Read();
    std::vector<std::string> gene_layout = GetGeneLayout(input_record.seq, delimiter);
    
    while (!(input_record = input_file->Read()).Empty()) {
        auto co_occurrences = GetCoOccurrences(gene_layout, input_record.seq, delimiter);

        if (co_occurrences.second.empty()) {
            // This gene doesn't interact with any other gene. Skipping it.
            if (arguments.verbose_output)
                fprintf(stdout, "The input table does not contain interactions for gene '%s'\n", co_occurrences.first.c_str());
            
            continue;
        }
        
        SequenceRecord output_record;
        for (const auto& interacting_gene: co_occurrences.second) {
            output_record.seq = co_occurrences.first;
            output_record.seq += delimiter;
            output_record.seq += interacting_gene;
            out_file->Write(output_record);
        }
    }
    fprintf(stdout, "The output file is located at \'%s\'\n", arguments.output_file_path.c_str());
    return 0;
}
