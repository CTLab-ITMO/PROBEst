# MIT License
#
# Copyright (c) 2025 CTLab-ITMO
#
# Authors: Daniil Smutin, Aleksandr Serdiukov, Vitalii Dravgelis, Artem Ivanov,
# Aleksei Zabashta, Sergey Muravyov, and the CTLab-ITMO university team.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from jsonschema import validate, ValidationError
import json
import argparse
import sys
import os
import pathlib

schema = {
  "type": "object",
  "properties": {
    "article_data": {
      "type": "object",
      "properties": {
        "probe_samples": {
          "type": "object",
          "properties": {
            "doi": { "type": "string" },
            "sample_group": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "probe_group": {
                    "type": "object",
                    "properties": {
                      "id_group": { "type": "integer" },
                      "experiment_ids": {
                        "type": "array",
                        "items": { "type": "integer" }
                      },
                      "results": {
                        "type": "object",
                        "properties": {
                          "result": {
                            "type": "array",
                            "items": {
                              "type": "object",
                              "properties": {
                                "id_exp": { "type": "integer" },
                                "outcome": { "type": "boolean" },
                                "fluorescence": { "type": "number" }
                              },
                              "required": ["id_exp", "outcome", "fluorescence"]
                            }
                          }
                        },
                        "required": ["result"]
                      },
                      "sequences": {
                        "type": "object",
                        "properties": {
                          "related_sequences": {
                            "type": "object",
                            "properties": {
                              "related_sequence": {
                                "type": "array",
                                "items": { "type": "string" }
                              }
                            },
                            "required": ["related_sequence"]
                          },
                          "probe_sequences": {
                            "type": "object",
                            "properties": {
                              "probes": {
                                "type": "object",
                                "properties": {
                                  "probe": {
                                    "type": "array",
                                    "items": {
                                      "type": "object",
                                      "properties": {
                                        "id_probe": { "type": "integer" },
                                        "probe_sequence": { "type": "string", "pattern": "^[ATGCatgc]{6,}$" },
                                        "modification": {
                                          "type": "array",
                                          "items": {
                                            "type": "object",
                                            "properties": {
                                              "modification_pos": { "type": "integer" },
                                              "modification_type": { "type": "string" }
                                            },
                                            "required": ["modification_pos", "modification_type"]
                                          }
                                        }
                                      },
                                      "required": ["id_probe", "probe_sequence", "modification"]
                                    }
                                  }
                                },
                                "required": ["probe"]
                              }
                            },
                            "required": ["probes"]
                          }
                        },
                        "required": ["related_sequences", "probe_sequences"]
                      }
                    },
                    "required": ["id_group", "experiment_ids", "results", "sequences"]
                  }
                },
                "required": ["probe_group"]
              }
            }
          },
          "required": ["doi", "sample_group"]
        },
        "hybridization_experiments": {
          "type": "object",
          "properties": {
            "doi": { "type": "string" },
            "experiment": {
              "type": "array",
              "items": {
                "type": "object",
                "properties": {
                  "id_exp": { "type": "integer" },
                  "raw_description": { "type": "string" },
                  "type": { "type": "string" },
                  "organism": { "type": "string" },
                  "rna_impurities": { "type": "boolean" },
                  "annealing": { "type": "boolean" },
                  "ph": { "type": "number" },
                  "concentrations": {
                    "type": "object",
                    "properties": {
                      "concentration": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "dna_rna_concentration": { "type": "string" },
                            "id_probe": { "type": "integer" },
                            "concentration_article": { "type": "string" },
                            "concentration_si": { "type": "number" }
                          },
                          "required": ["dna_rna_concentration", "id_probe", "concentration_article", "concentration_si"]
                        }
                      }
                    },
                    "required": ["concentration"]
                  },
                  "string_parameters": {
                    "type": "object",
                    "properties": {
                      "parameter": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "temperature": { "type": "string" },
                            "tris": { "type": "string" },
                            "na": { "type": "string" },
                            "k": { "type": "string" },
                            "mg": { "type": "string" },
                            "dmso": { "type": "string" }
                          },
                          "required": ["temperature", "tris", "na", "k", "mg", "dmso"]
                        }
                      }
                    },
                    "required": ["parameter"]
                  },
                  "numeric_parameters": {
                    "type": "object",
                    "properties": {
                      "numeric_parameter": {
                        "type": "array",
                        "items": {
                          "type": "object",
                          "properties": {
                            "temperature": { "type": "number" },
                            "tris": { "type": "number" },
                            "na": { "type": "number" },
                            "k": { "type": "number" },
                            "mg": { "type": "number" },
                            "dmso": { "type": "number" }
                          },
                          "required": ["temperature", "tris", "na", "k", "mg", "dmso"]
                        }
                      }
                    },
                    "required": ["numeric_parameter"]
                  }
                },
                "required": ["id_exp", "raw_description", "type", "organism", "rna_impurities", "annealing", "ph", "concentrations", "string_parameters", "numeric_parameters"]
              }
            }
          },
          "required": ["doi", "experiment"]
        }
        },
          "required": ["probe_samples", "hybridization_experiments"]
    },
  },
  "required": ["article_data"]
}



def validate_json(json_data, schema):
    try:
        validate(instance=json_data, schema=schema)
        return True
    except ValidationError as e:
        print(f"Validation error: {e.message}")
        return False

def load_json_file(file_path):
    with open(file_path, "r") as file:
        return json.load(file)


def main(cmdline=None):
    if cmdline is None:
        cmdline=sys.argv[1:]
    parser = argparse.ArgumentParser(
        description='This scripts validate JSON from LLM against determined schema with data types and RegExp')
    parser.add_argument('file_path', type=str, help=
                        "Path to JSON file to validate")

    args = parser.parse_args(cmdline)
    file_path = args.file_path
    if not os.path.isabs(file_path):
            file_path = pathlib.Path(os.getcwd()) / pathlib.Path(file_path)
    result = validate_json(load_json_file(file_path), schema)
    if result:
        print('Pass')

if __name__ == '__main__':
    main()