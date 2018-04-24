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


#ifndef CoOccurrenceReporter_hpp
#define CoOccurrenceReporter_hpp

#include "../libgene/source/utils/Tokenizer.hpp"

#include <string>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>

//
// Returns a vector of individual genes symbols.
//
// 'raw_line' is the first line of the input file which contains the column
// information. The first column is dropped since it doesn't contain any useful
// information.
//
std::vector<std::string> GetGeneLayout(const std::string& raw_line,
                                       char delimiter)
{
    std::vector<std::string> gene_layout;
    Tokenizer splitter(raw_line, delimiter);
    
    splitter.ReadNext();
    std::string id_column = splitter.GetNextToken();
    
    while (splitter.ReadNext()) {
        std::string gene = splitter.GetNextToken();
        gene_layout.emplace_back(std::move(gene));
    }
    return gene_layout;
}

//
// Returns a pair of gene symbol and its interactions with other genes. The
// order is determined by decreasing absolute value of each one's strength.
//
auto GetCoOccurrences(const std::vector<std::string>& gene_layout,
                      const std::string& line,
                      const char delimiter)
{
    std::pair<std::string, std::vector<std::string>> result;
    
    Tokenizer splitter(delimiter);
    splitter.SetText(line);
    
    splitter.ReadNext();
    std::string gene = splitter.GetNextToken();
    result.first = std::move(gene);
    
    std::vector<std::pair<std::string, double>> cooccurrences;
    
    int index = 0;
    while (splitter.ReadNext()) {
        std::string value = splitter.GetNextToken();
        if (!value.empty()) {
            // There is some interaction between this gene (in schema) and the
            // one we're currently at
            cooccurrences.push_back(std::make_pair(gene_layout[index],
                                                   std::abs(std::stod(value))));
        }
        ++index;
    }

    // Sort interactions based on their strength
    std::sort(cooccurrences.begin(), cooccurrences.end(),
              [](auto&& arg1, auto&& arg2) {
                  return arg1.second > arg2.second;
              });
    
    result.second.reserve(cooccurrences.size());
    // Write interactions out in the correct order, dropping the strength
    for (auto& cooccurrence: cooccurrences) {
        result.second.emplace_back(std::move(cooccurrence.first));
    }
    return result;
}

#endif /* CoOccurrenceReporter_hpp */
