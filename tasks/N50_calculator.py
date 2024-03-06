#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

length_file = sys.argv[1]

def calculate_N50(list_of_lengths):
    """Calculate N50 for a sequence of numbers.

    Args:
        list_of_lengths (list): List of numbers.

    Returns:
        float: N50 value.

    """
    tmp = []
    for tmp_number in set(list_of_lengths):
        tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()

    if (len(tmp) % 2) == 0:
        median = (tmp[int(len(tmp) / 2) - 1] + tmp[int(len(tmp) / 2)]) / 2
    else:
        median = tmp[int(len(tmp) / 2)]

    return median


def calculate_N90(list_of_lengths):
    """Calculate N50 for a sequence of numbers.

    Args:
        list_of_lengths (list): List of numbers.

    Returns:
        float: N50 value.

    """
    tmp = []
    for tmp_number in set(list_of_lengths):
        tmp += [tmp_number] * list_of_lengths.count(tmp_number) * tmp_number
    tmp.sort()

    if (len(tmp) % (1/.1)) == 0:
        median = (tmp[int(len(tmp) / (1/.1)) - 1] + tmp[int(len(tmp) / (1/.1))]) / 2
    else:
        median = tmp[int(len(tmp) / (1/.1))]

    return median





length_list = open(length_file, "r").read().splitlines()
length_list = [int(i) for i in length_list]

print("N50\tN90")
print(str(calculate_N50(length_list)) + "\t" + str(calculate_N90(length_list)))

