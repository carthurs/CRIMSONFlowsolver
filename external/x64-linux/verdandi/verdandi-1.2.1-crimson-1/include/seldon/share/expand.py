# Copyright (C) 2001-2009 Vivien Mallet
#
# This file is part of the linear-algebra library Seldon,
# http://seldon.sourceforge.net/.
#
# Seldon is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Seldon. If not, see http://www.gnu.org/licenses/.


import sys, os, re


#######################
# PARSES COMMAND LINE #
#######################


if len(sys.argv) != 3:
    print "Please provide exactly two input files."
    sys.exit(1)

template_definition = sys.argv[1]
input_file = sys.argv[2]
if len(input_file) < 9 or input_file[:9] != "template-":
    print "Input filename should start with \"template-\"."
    sys.exit(1)

if not os.path.isfile(input_file):
    print "\"" + input_file + "\" is not a file!"
    sys.exit(1)

output_file = input_file[9:]
if output_file == "":
    print "Error: output filename is empty."
    sys.exit(1)


####################
# USEFUL FUNCTIONS #
####################


def combine(list_of_list):
    """
    Generates all combinations (x, y, ...) where x is an element of
    list_of_list[0], y an element of list_of_list[1], ...

    \param list_of_list the list of lists of possible values.
    \return The list of all possible combinations (therefore, a list of
    lists).
    """
    if len(list_of_list) == 0:
        return []
    elif len(list_of_list) == 1:
        return [[x] for x in list_of_list[0]]
    elif len(list_of_list) == 2:
        return [[x, y] for x in list_of_list[0] for y in list_of_list[1]]
    else:
        middle = len(list_of_list) / 2
        return [x + y for x in combine(list_of_list[:middle])
                for y in combine(list_of_list[middle:])]


def expand_string(input_string, replacement_list, recursive = True):
    """
    Expands markups in a string, based on all possible values that may take
    these markups.
    
    \param input_string the input string in which the markups are to be
    replaced.
    \param replacement_list a dictionary which markups in keys and a list of
    possible values for each key.
    \param recursive are recursive expansions enabled?
    \return The list of all possible strings.
    """
    key_list = replacement_list.keys()
    # key_list should be sorted from the longest string to the smallest
    # string, so that the replacements do not interfere with each other.
    key_list.sort(key = len)
    key_list.reverse()
    # Active markups.
    active = []
    tmp_string = input_string
    for markup in key_list:
        active.append(markup in tmp_string)
        # If any active markup is not to be expanded, an empty list should be
        # returned.
        if markup in tmp_string and replacement_list[markup] == []:
            return []
        tmp_string = tmp_string.replace(markup, "")
    # Is any replacement needed?
    if any(active):
        output_list = []
        # Selects the replacements actually needed.
        active_key = [key for key, act in zip(key_list, active) if act]
        active_replacement = [replacement_list[key] for key in active_key]
        # Generates all possible combinations.
        active_replacement = combine(active_replacement)
        for replacement in active_replacement:
            output_string = input_string
            # Expands the markups one by one.
            for key, value in zip(active_key, replacement):
                output_string = output_string.replace(key, value)
            output_list.append(output_string)
        # Recursion may be needed to complete the expansion.
        if recursive:
            raw_list = output_list
            output_list = []
            for item in raw_list:
                output_list += expand_string(item, replacement_list, True)
        return output_list
    else:
        return [input_string]


def read_markup_list(filename):
    """
    Extracts a replacement list from a file.

    \param filename name of the file that defines the markups.
    \return A dictionary with the markups and their possible values.
    """
    # File contents.
    file_stream = open(filename)
    file_line = file_stream.readlines()
    file_stream.close()
    Nline = len(file_line)

    # Searches for list definition: (define ...).
    define_list = []
    i = 0
    while i != Nline:
        if "(define" in file_line[i].strip():
            line = ""
            while True:
                line += file_line[i].strip() + " "
                if ")" in line:
                    define_list.append(line.strip())
                    i += 1
                    break
                elif i == Nline:
                    raise Exception, "Syntax error: \"define\" not closed."
                i += 1
                if i == Nline:
                    break
        else:
            i += 1

    # Now parses the definitions.
    dictionary = {}
    for line in define_list:
        if len(line) < 13 or line[:8] != "(define " or line[-2:] != ");":
            raise Exception, "Syntax error in \"define\":\n" + line
        if "(" in line[8:-2]:
            raise Exception, "Syntax error in \"define\":\n" + line
        line_list = re.split("[ \n\t,:=]+", line[8:-2].strip())
        if line_list != [] and line_list[-1] == "":
            # Removes the empty string.
            line_list = line_list[:-1]
        if line_list == []:
            raise Exception, "Syntax error in \"define\":\n" + line
        if len(line_list) == 1:
            # Empty: no replacement.
            dictionary[line_list[0]] = []
        else:
            dictionary[line_list[0]] = line_list[1:]

    return dictionary
        

def read_expand(filename, replacement_list):
    """
    Applies a replacement list to a file. Instructions "(define ...)" are
    filtered out.

    \param filename name of the file to proceed with.
    \return A list of the line with the markups expanded.
    """
    # File contents.
    file_stream = open(filename)
    file_line = file_stream.readlines()
    file_stream.close()
    Nline = len(file_line)

    # Searches for code lines, that is, all lines but those with a definition:
    # (define ...).
    code_line = []
    i = 0
    while i != Nline:
        if "(define" in file_line[i].strip():
            while True:
                if ")" in file_line[i]:
                    i += 1
                    break
                elif i == Nline:
                    raise Exception, "Syntax error: \"define\" not closed."
                i += 1
                if i == Nline:
                    break
        else:
            code_line.append(file_line[i])
            i += 1

    # Now proceeds with the replacements.
    output_line = []
    for line in code_line:
        output_line += expand_string(line, replacement_list)

    return output_line
        

###################
# FILE GENERATION #
###################


replacement_list = read_markup_list(template_definition)

output_stream = open(output_file, "w")
output_stream.write("".join(read_expand(input_file, replacement_list)))
output_stream.close()
