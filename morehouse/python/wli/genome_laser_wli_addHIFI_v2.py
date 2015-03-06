#!/usr/bin/python
#######################################################################################
# Guoxing Fu Jan 28, 2015
# 
#######################################################################################

import os, glob, subprocess, random, operator, time, sys, copy, shutil, fnmatch
from optparse import OptionParser


class parameters:
    def __init__(self):
        self.person_dict = {}
        self.rsID_dict = {}
        self.pos_list = {}
        self.father_list = []
        self.mather_list = []
        self.children_list = []
        self.children_dict = {}

        self.fragment_dict = {}
        self.fragment_list = []

        self.children_hap_file = "child_hap.txt"


class persons:
    def __init__(self):
        self.ID = ""
        self.father = ""
        self.mather = ""
        self.children = {}
        self.genotype_dict = {}
        self.haplotype = {}
        self.hetero_pos_list = []


class fragments:
    def __init__(self):
        self.ID = ""
        self.start = ""
        self.end = ""
        self.length = ""


def list_to_line(list):
    line = ""
    for a in list:
        line += str(a).strip() + "\t"
    return line.strip()


def sort_dict_by_key(input_dict):
    sorted_list = []
    sorted_list = [x for x in input_dict.iteritems()]
    sorted_list.sort(key=lambda x: x[0])
    return sorted_list


def load_pedi(pedi_name):
    with open(pedi_name, "r") as pedi_file:
        for line in pedi_file:
            if line.startswith("PersonID"):
                title_info = list_to_line(line.strip().split())
            else:
                try:
                    elements = line.strip().split()
                    person = persons()
                    person.ID = elements[0]
                    person.father = elements[1]
                    person.mather = elements[2]
                    if person.father != "N/A" and person.father in parameter.person_dict and person.ID not in \
                            parameter.person_dict[person.father].children:
                        parameter.person_dict[person.father].children[person.ID] = ""
                    if person.mather != "N/A" and person.mather in parameter.person_dict and person.ID not in \
                            parameter.person_dict[person.mather].children:
                        parameter.person_dict[person.mather].children[person.ID] = ""
                    if person.ID not in person_dict:
                        parameter.person_dict[person.ID] = person
                except:
                    print "error in ", line, pedi_name


def load_geno(geno_name):
    with open(geno_name, "r") as geno_file:
        for line in geno_file:
            if line.startswith("rs#"):
                ID_list = line.strip().split()[2:]
            else:
                try:
                    elements = line.strip().split()
                    rsID = elements[0]
                    position = int(elements[1])
                    parameter.rsID_dict[position] = rsID
                    genotype = elements[2:]
                    for index, ID in enumerate(ID_list):
                        parameter.person_dict[ID].genotype_dict[position] = genotype[index]
                except:
                    print "error in geno", line, pedi_name
                    # pass


def parents_to_children():
    for ID in parameter.person_dict:
        person = parameter.person_dict[ID]
        if person.father != "N/A" and person.mather != "N/A":

            pos_list = person.genotype_dict.keys()
            pos_list.sort()
            for pos in pos_list:
                c_geno = person.genotype_dict[pos]
                f_geno = parameter.person_dict[person.father].genotype_dict[pos]
                m_geno = parameter.person_dict[person.mather].genotype_dict[pos]

                if True:
                    f_set = set(f_geno)
                    m_set = set(m_geno)
                    c_set = set(c_geno)

                    if len(c_set) == 1 and c_geno != "NN":
                        if (c_geno[0] == f_geno[0] or c_geno[0] == f_geno[1]) and (
                                c_geno[0] == m_geno[0] or c_geno[0] == m_geno[1]):
                            parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
                        """
                        # keep this part. for non-NN c_geno, has discrepancy with f_geno or m_geno. uesful for future.
                        else:
                            if f_geno == "NN" and m_geno == "NN":
                                parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
                            elif f_geno == "NN" and m_geno != "NN" and (c_geno[0] == m_geno[0] or c_geno[0] == m_geno[1]):
                                parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
                            elif f_geno != "NN" and m_geno == "NN" and (c_geno[0] == f_geno[0] or c_geno[0] == f_geno[1]):
                                parameter.person_dict[ID].haplotype[pos] = (c_geno[0], c_geno[1])
                            else:
                                #pass
                                print "child geno not in parent geno", ID, pos, f_geno, m_geno, c_geno
                                #sys.exit(1)
                        """
                    else:
                        if f_geno != "NN":
                            if m_geno != "NN":
                                if c_geno != "NN":
                                    if len(f_set) == 1:
                                        if len(m_set) == 1:
                                            parameter.person_dict[ID].haplotype[pos] = (f_geno[0], m_geno[1])
                                        elif len(m_set) != 1:
                                            cf_hap = f_geno[0]
                                            cm_hap = c_geno[1] if c_geno[0] == f_geno[0] else c_geno[0]
                                            parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
                                    elif len(f_set) != 1:
                                        if len(m_set) == 1:
                                            cf_hap = c_geno[1] if c_geno[0] == m_geno[0] else c_geno[0]
                                            cm_hap = m_geno[0]
                                            parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
                                        elif len(m_set) != 1:
                                            # all hetero, cannot determine
                                            parameter.person_dict[ID].haplotype[pos] = ("X", "X")

                                elif c_geno == "NN":
                                    if len(f_set) == 1:
                                        if len(m_set) == 1:
                                            parameter.person_dict[ID].haplotype[pos] = (f_geno[0], m_geno[1])
                                        elif len(m_set) != 1:
                                            parameter.person_dict[ID].haplotype[pos] = (f_geno[0], "N")
                                    elif len(f_set) != 1:
                                        if len(m_set) == 1:
                                            parameter.person_dict[ID].haplotype[pos] = ("N", m_geno[0])
                                        elif len(m_set) != 1:
                                            parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                                else:
                                    print "error 1", pos

                            elif m_geno == "NN":
                                if c_geno != "NN":
                                    if len(f_set) == 1:
                                        cf_hap = f_geno[0]
                                        cm_hap = c_geno[1] if c_geno[0] == f_geno[0] else c_geno[0]
                                        parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
                                    elif len(f_set) != 1:
                                        parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                                elif c_geno == "NN":
                                    if len(f_set) == 1:
                                        parameter.person_dict[ID].haplotype[pos] = (f_geno[0], "N")
                                    elif len(f_set) != 1:
                                        parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                            else:
                                print "error 2", pos

                        elif f_geno == "NN":
                            if m_geno != "NN":
                                if c_geno != "NN":
                                    if len(m_set) == 1:
                                        cf_hap = c_geno[1] if c_geno[0] == m_geno[0] else c_geno[0]
                                        cm_hap = m_geno[0]
                                        parameter.person_dict[ID].haplotype[pos] = (cf_hap, cm_hap)
                                    elif len(m_set) != 1:
                                        parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                                if c_geno == "NN":
                                    if len(m_set) == 1:
                                        parameter.person_dict[ID].haplotype[pos] = ("N", m_geno[0])
                                    elif len(m_set) != 1:
                                        parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                            if m_geno == "NN":
                                parameter.person_dict[ID].haplotype[pos] = ("N", "N")
                        else:
                            print "3", pos
                    try:
                        # if len(c_set) != 1:
                        #print f_geno, m_geno, c_geno, parameter.person_dict[ID].haplotype[pos][0], parameter.person_dict[ID].haplotype[pos][1]
                        pass
                    except:
                        pass
                        # print "error 4", pos, ID, f_geno, m_geno, c_geno
                    #print "child geno not in parent geno", ID, pos, f_geno, m_geno, c_geno
                    #sys.exit(1)


def output_child_hap():
    pos_list = parameter.rsID_dict.keys()
    pos_list.sort()

    with open(parameter.children_hap_file, "w") as c_hap_file:
        print >> c_hap_file, "rs_ID", "pos",
        for id in parameter.children_list:
            print >> c_hap_file, id + "_F", id + "_M",
        print >> c_hap_file, ""

        for pos in pos_list:
            print >> c_hap_file, parameter.rsID_dict[pos], pos,
            for child_id in parameter.children_list:
                print >> c_hap_file, parameter.person_dict[child_id].haplotype[pos][0], \
                parameter.person_dict[child_id].haplotype[pos][1],
            print >> c_hap_file, ""


def prepare_id_list():
    # children list
    for ID in parameter.person_dict.keys():
        person = parameter.person_dict[ID]
        if len(person.children) > 0:
            parameter.children_list.extend(list(person.children.keys()))
    # print parameter.children_list

    parameter.children_list = list(set(parameter.children_list))
    parameter.children_list.sort()
    # print parameter.children_list

    # parent list
    for ID in parameter.person_dict.keys():
        person = parameter.person_dict[ID]
        if person.father != "N/A":
            parameter.father_list.append(person.father)
        if person.mather != "N/A":
            parameter.mather_list.append(person.mather)

    parameter.father_list = list(set(parameter.father_list))
    parameter.father_list.sort()
    #print parameter.father_list

    parameter.mather_list = list(set(parameter.mather_list))
    parameter.mather_list.sort()
    #print parameter.mather_list

    #parameter.pos_list = [int(x) for x in parameter.rsID_dict.keys()]
    parameter.pos_list = parameter.rsID_dict.keys()
    parameter.pos_list.sort()


def sort_fragment(f_id):
    # f_id = "1NA19702"

    children_list = parameter.person_dict[f_id].children.keys()
    children_list.sort()
    #print children_list

    fragment_startpos_dict = {}

    total_fragment = 0
    unique_fragment_total = 0

    for child_id in children_list:
        unique_fragment = []
        list = parameter.fragment_dict[child_id]
        for fragment in list:
            total_fragment += 1
            in_unique_list = False
            #print fragment.ID, fragment.start, fragment.end
            start_pos = fragment.start
            end_pos = fragment.end
            for uni_frangment in unique_fragment:
                if uni_frangment.start_pos == start_pos and uni_frangment.end_pos == end_pos:
                    in_unique_list = True
                    break

            if not in_unique_list:
                unique_fragment_total += 1
                if start_pos not in fragment_startpos_dict:
                    fragment_startpos_dict[int(start_pos)] = []
                fragment_startpos_dict[int(start_pos)].append(fragment)

    #print "total_fragment", total_fragment
    #print "unique_fragment_total", unique_fragment_total

    fragment_startpos_sorted_list = sort_dict_by_key(fragment_startpos_dict)
    """

    for child_id in children_list:
        list = parameter.fragment_dict[child_id]
        for fragment in list:
            print fragment.ID, fragment.start, fragment.end
            start_pos = fragment.start
            end_pos = fragment.end
            if start_pos not in fragment_startpos_dict:
                fragment_startpos_dict[int(start_pos)] = []
            fragment_startpos_dict[int(start_pos)].append(fragment)

    fragment_startpos_sorted_list = sort_dict_by_key(fragment_startpos_dict)


    for data in fragment_startpos_sorted_list:
        start_pos = data[0]
        fragment_list = data[1]

        print start_pos,
        for fragment in fragment_list:
            print fragment.ID, fragment.start
    """
    return fragment_startpos_sorted_list


def children_to_parents_max_fragment():
    f_id = "1NA19702"

    children_list = parameter.person_dict[f_id].children.keys()
    children_list.sort()
    # print children_list

    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[f_id].genotype_dict[pos]
        if f_geno != "NN" and f_geno[0] != f_geno[1]:
            parameter.person_dict[f_id].hetero_pos_list.append(pos)

    #print "hetero pos_list size", len(parameter.person_dict[f_id].hetero_pos_list)

    #for child_id in children_list:

    compare_child_hap(children_list[0], children_list[1])
    compare_child_hap(children_list[0], children_list[2])
    compare_child_hap(children_list[1], children_list[2])

    max_length_fragment = ""
    max_length = 0

    for child_id in children_list:
        list = parameter.fragment_dict[child_id]
        for fragment in list:
            if max_length == 0:
                max_lsame_to_Aength = fragment.length
                max_length_fragment = fragment
            else:
                if fragment.length > max_length:
                    max_length = fragment.length
                    max_length_fragment = fragment
                #print fragment.ID, fragment.start, fragment.end, fragment.length
    #print max_length, max_length_fragment.ID

    temp_parent_hap_A = {}
    temp_parent_hap_B = {}
    temp_parent_hap_list = [temp_parent_hap_A, temp_parent_hap_B]

    for pos in parameter.pos_list:
        if pos >= max_length_fragment.start and pos <= max_length_fragment.end:
            temp_parent_hap_A[pos] = parameter.person_dict[max_length_fragment.ID].haplotype[pos][0]
        #print pos, temp_parent_hap_A[pos]
    #parameter.fragment_dict[child_id].remove(max_length_fragment)

    child_id = children_list[1]
    #if True:
    for child_id in children_list:
        if True:
            #if child_id != max_length_fragment.ID:
            list = parameter.fragment_dict[child_id]
            fragment = list[0]
            if True:
                #while len(parameter.fragment_dict[child_id]) > 0:
                #print "parameter.fragment_dict[child_id] size", len(parameter.fragment_dict[child_id])
                for fragment in list:
                    #current_start, current_end = get_parent_hap_end(temp_parent_hap_A)
                    status, temp_parent_hap_A = update_parent_hap(fragment, temp_parent_hap_A)
                    if status == "same_to_A":
                        parameter.fragment_dict[child_id].remove(fragment)
                    #print fragment.ID, fragment.start, "removed", len(parameter.fragment_dict[child_id])

                    elif status == "same_to_B":
                        status, temp_parent_hap_B = update_parent_hap(fragment, temp_parent_hap_B)
                    elif status == "mix":
                        pass
    #print "temp_parent_hap_A", len(temp_parent_hap_A)
    #print "temp_parent_hap_B", len(temp_parent_hap_B)


    count = 0
    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[f_id].genotype_dict[pos]
        if f_geno[0] == f_geno[1]:
            if pos not in temp_parent_hap_A:
                temp_parent_hap_A[pos] = f_geno[0]
            if pos not in temp_parent_hap_B:
                temp_parent_hap_B[pos] = f_geno[0]
        else:
            if pos in temp_parent_hap_A and pos not in temp_parent_hap_B:
                temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]

            elif pos not in temp_parent_hap_A and pos in temp_parent_hap_B:
                temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]

            elif pos not in temp_parent_hap_A and pos not in temp_parent_hap_B:
                temp_parent_hap_A[pos] = "N"
                temp_parent_hap_B[pos] = "N"
                count += 1

    #print "count", count
    """
    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[f_id].genotype_dict[pos]
        if f_geno[0] != f_geno[1]:
            print pos,
            if pos in temp_parent_hap_A:
                print temp_parent_hap_A[pos],
            else:
                print " ",
            if pos in temp_parent_hap_B:
                print temp_parent_hap_B[pos]
            else:
                print " "

    for pos in temp_parent_hap_A:
        #print pos, temp_parent_hap_A[pos],
        if pos not in temp_parent_hap_B:
            print pos, temp_parent_hap_A[pos]
    """

# print "temp_parent_hap_A final", len(temp_parent_hap_A)
# print "temp_parent_hap_B final", len(temp_parent_hap_B)

def children_to_parents(p_id, p_code):

    children_list = parameter.person_dict[p_id].children.keys()
    children_list.sort()
    #print children_list

    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[p_id].genotype_dict[pos]
        if f_geno != "NN" and f_geno[0] != f_geno[1]:
            parameter.person_dict[p_id].hetero_pos_list.append(pos)

    #print "hetero pos_list size", len(parameter.person_dict[f_id].hetero_pos_list)

    #for child_id in children_list:

    compare_child_hap(p_id, p_code, children_list[0], children_list[1])
    compare_child_hap(p_id, p_code, children_list[0], children_list[2])
    compare_child_hap(p_id, p_code, children_list[1], children_list[2])

    temp_parent_hap_A = {}
    temp_parent_hap_B = {}

    fragment_startpos_sorted_list = sort_fragment(p_id)
    status = ""
    for data in fragment_startpos_sorted_list:
        start_pos = data[0]
        fragment_list = data[1]

        #print start_pos
        for fragment in fragment_list:
            #print fragment.ID, fragment.start, fragment.end
            status, temp_parent_hap_A = update_parent_hap(p_id, fragment, temp_parent_hap_A, p_code)
            if status == "same_to_A":
                pass
            elif status == "same_to_B":
                status, temp_parent_hap_B = update_parent_hap(p_id, fragment, temp_parent_hap_B, p_code)
            elif status == "mix":
                pass

    #print "temp_parent_hap_A", len(temp_parent_hap_A)
    #print "temp_parent_hap_B", len(temp_parent_hap_B)

    count = 0
    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[p_id].genotype_dict[pos]
        if f_geno[0] == f_geno[1]:
            #if pos not in temp_parent_hap_A:
            temp_parent_hap_A[pos] = f_geno[0]
            #if pos not in temp_parent_hap_B:
            temp_parent_hap_B[pos] = f_geno[0]
        else:
            if pos in temp_parent_hap_A and pos not in temp_parent_hap_B:
                temp_parent_hap_B[pos] = f_geno[0] if temp_parent_hap_A[pos] == f_geno[1] else f_geno[1]

            elif pos not in temp_parent_hap_A and pos in temp_parent_hap_B:
                temp_parent_hap_A[pos] = f_geno[0] if temp_parent_hap_B[pos] == f_geno[1] else f_geno[1]

            elif pos not in temp_parent_hap_A and pos not in temp_parent_hap_B:
                temp_parent_hap_A[pos] = "N"
                temp_parent_hap_B[pos] = "N"
                count += 1

    #print "count N", count


    for pos in parameter.pos_list:
        if pos not in parameter.person_dict[p_id].haplotype:
            parameter.person_dict[p_id].haplotype[pos] = ["", ""]

        if pos in temp_parent_hap_A:
            parameter.person_dict[p_id].haplotype[pos][0] = temp_parent_hap_A[pos]
        else:
            parameter.person_dict[p_id].haplotype[pos][0] = "N"
        if pos in temp_parent_hap_B:
            parameter.person_dict[p_id].haplotype[pos][1] = temp_parent_hap_B[pos]
        else:
            parameter.person_dict[p_id].haplotype[pos][1] = "N"

	"""
    # to output parent hap data
    with open(p_id + "_hap.txt", "w") as x_file:
        for pos in parameter.pos_list:
            f_geno = parameter.person_dict[p_id].genotype_dict[pos]
            #if f_geno[0] != f_geno[1]:
            if True:
                print >> x_file, pos,
                if pos in temp_parent_hap_A:
                    print >> x_file, temp_parent_hap_A[pos],
                else:
                    print >> x_file, "N",
                if pos in temp_parent_hap_B:
                    print >> x_file, temp_parent_hap_B[pos]
                else:
                    print >> x_file, "N"

                #print "temp_parent_hap_A final", len(temp_parent_hap_A)
                #print "temp_parent_hap_B final", len(temp_parent_hap_B)
    """


def output_parent_hap():
    parent_id_list = []
    parent_id_list.extend(parameter.father_list)
    parent_id_list.extend(parameter.mather_list)
    parent_id_list.sort()

    with open("parent_hap.txt", "w") as parent_hap:
        print >> parent_hap, "rs#", "pos",
        for id in parent_id_list:
            print >> parent_hap, parameter.person_dict[id].ID + "_A", parameter.person_dict[id].ID + "_B",
        print >> parent_hap, ""

        for pos in parameter.pos_list:
            print >> parent_hap, parameter.rsID_dict[pos], pos,
            for id in parent_id_list:
                print >> parent_hap, parameter.person_dict[id].haplotype[pos][0], \
                parameter.person_dict[id].haplotype[pos][1],
            print >> parent_hap, ""


def update_parent_hap(f_id, fragment, temp_parent_hap_A, p_index):
    if len(temp_parent_hap_A) == 0:
        # for temp_parent_hap_B
        for pos in parameter.pos_list:
            if pos >= fragment.start and pos <= fragment.end:
                temp_parent_hap_A[pos] = parameter.person_dict[fragment.ID].haplotype[pos][p_index]
        return "same_to_A", temp_parent_hap_A
    else:
        current_start, current_end = get_parent_hap_end(temp_parent_hap_A)

        overlap_start = 0
        overlap_end = 0

        if fragment.end <= current_start or fragment.start >= current_end:
            return "no_overlap", temp_parent_hap_A
        else:
            if fragment.start < current_start and fragment.end < current_end:
                overlap_start = current_start
                overlap_end = fragment.end
            elif fragment.start > current_start and fragment.end < current_end:
                overlap_start = fragment.start
                overlap_end = fragment.end
            else:
                overlap_start = fragment.start
                overlap_end = current_end

        #print "current_start, current_end, overlap_start, overlap_end", current_start, current_end, overlap_start, overlap_end

        same = 0
        not_same = 0
        #f_id = "1NA19702"
        for pos in parameter.person_dict[f_id].hetero_pos_list:
            if pos >= overlap_start and pos <= overlap_end:
                if parameter.person_dict[fragment.ID].haplotype[pos][p_index] != "X" and temp_parent_hap_A[pos] != "X" \
                        and parameter.person_dict[fragment.ID].haplotype[pos][p_index] != "N" and temp_parent_hap_A[
                    pos] != "N":
                    #print "xxxx", pos, parameter.person_dict[fragment.ID].haplotype[pos][0], temp_parent_hap[pos]
                    if parameter.person_dict[fragment.ID].haplotype[pos][p_index] == temp_parent_hap_A[pos]:
                        same += 1
                    else:
                        not_same += 1

        #print "before", len(temp_parent_hap_A)
        percentage = float(same) / (same + not_same + 1)
        #print percentage
        if percentage > 0.9:
            for pos in parameter.person_dict[f_id].hetero_pos_list:
                if pos >= fragment.start and pos <= fragment.end and pos not in temp_parent_hap_A:
                    temp_parent_hap_A[pos] = parameter.person_dict[fragment.ID].haplotype[pos][p_index]
            #print "added snp", len(temp_parent_hap_A)
            return "same_to_A", temp_parent_hap_A
        elif percentage < 0.1:
            return "same_to_B", temp_parent_hap_A
        else:
            #print fragment.ID, fragment.start, fragment.end, fragment.length, percentage
            return "mix", temp_parent_hap_A


def get_parent_hap_end(temp_parent_hap):
    pos_list = temp_parent_hap.keys()
    pos_list.sort()
    return (pos_list[0], pos_list[-1])


def compare_child_hap(p_id, p_code, child_ID_1, child_ID_2):
    child_fragment_dict = {}
    child_fragment_dict[child_ID_1] = []
    child_fragment_dict[child_ID_2] = []

    fragment_list = []

    break_point_dict = {}
    break_point_list = []

    child_1 = parameter.person_dict[child_ID_1]
    child_2 = parameter.person_dict[child_ID_2]

    child_fragment_dict[child_ID_1].append(parameter.pos_list[0])
    child_fragment_dict[child_ID_2].append(parameter.pos_list[0])

    compare = ""

    for pos in parameter.pos_list:
        f_geno = parameter.person_dict[p_id].genotype_dict[pos]

        if f_geno != "NN" and f_geno[0] != f_geno[1]:

            cf_hap_1 = child_1.haplotype[pos][p_code]
            cf_hap_2 = child_2.haplotype[pos][p_code]

            parameter.person_dict[p_id].hetero_pos_list.append(pos)

            if cf_hap_1 != "X" and cf_hap_2 != "X" and cf_hap_1 != "N" and cf_hap_2 != "N":

                if compare == "":
                    #print f_geno, cf_hap_1, cf_hap_2
                    compare = "same" if cf_hap_1 == cf_hap_2 else "different"
                #print "initial compare status:", compare

                #print "all", pos, cf_hap_1, cf_hap_2
                elif compare == "same":
                    if cf_hap_1 != cf_hap_2:
                        child_fragment_dict[child_ID_1].append(pos)
                        child_fragment_dict[child_ID_2].append(pos)
                        compare = "different"
                    #print "different", pos, cf_hap_1, cf_hap_2
                elif compare == "different":
                    if cf_hap_1 == cf_hap_2:
                        child_fragment_dict[child_ID_1].append(pos)
                        child_fragment_dict[child_ID_2].append(pos)
                        compare = "same"
                    #print "same", pos, cf_hap_1, cf_hap_2

    child_fragment_dict[child_ID_1].append(parameter.pos_list[-1])
    child_fragment_dict[child_ID_2].append(parameter.pos_list[-1])

    #print child_fragment_dict[child_ID_1]
    #print child_fragment_dict[child_ID_2]

    for id in child_ID_1, child_ID_2:
        for index, pos in enumerate(child_fragment_dict[id]):
            if index < len(child_fragment_dict[id]) - 1:
                fragment = fragments()
                fragment.ID = id
                fragment.start = child_fragment_dict[id][index]
                fragment.end = child_fragment_dict[id][index + 1]
                fragment.length = parameter.pos_list.index(fragment.end) - parameter.pos_list.index(fragment.start)

                #print fragment.ID, fragment.start, fragment.end, fragment.length

                if id not in parameter.fragment_dict:
                    parameter.fragment_dict[id] = []
                parameter.fragment_dict[id].append(fragment)


def genome_laser(pedi_name, geno_name):
    #pedi_name = "pedi.txt"
    #geno_name = "geno.txt"

    load_pedi(pedi_name)
    load_geno(geno_name)

    prepare_id_list()
    # print parameter.person_dict.keys()
    parents_to_children()
    output_child_hap()

    f_code = 0
    m_code = 0

    for f_id in parameter.father_list:
        children_to_parents(f_id, f_code)

    for m_id in parameter.mather_list:
        children_to_parents(m_id, m_code)

    output_parent_hap()


def get_args():
    desc = ""
    usage = ""
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option("-p", "--pedi", type="string", dest="pedi_name", help="Input pedi name", default="null")
    parser.add_option("-g", "--geno", type="string", dest="geno_file", help="Input file name", default="null")
    parser.add_option("-r", "--reference", type="string", dest="ref_file", help="Input reference name", default="null")
    (options, args) = parser.parse_args()

    return options

####################
# #def for HIFI# #
####################


def load_raw_data(hap_file):
    title = ""
    data = {}
    with open(hap_file, "r") as fp:
        for line in fp:
            if line != "":
                elements = line.strip().split()
                if line.startswith("rs_ID") or line.startswith("rs#") or line.startswith("rsID"):
                    title = elements
                else:
                    # print elements[1]
                    data[int(elements[1])] = elements
    return title, data


def find_index(a, b):
    index_list = []
    for str_b in b:
        if a in str_b:
            # print b.index(str_b)
            index_list.append(b.index(str_b))
    return index_list


def remove_sample_from_ref(line, rm_list):  # Deleting multiple elements from a list by index
    # new_line = [i for j, i in enumerate(line) if i not in rm_list]
    new_line = line
    for index in sorted(rm_list, reverse=True):
        del new_line[index]
    return new_line


def make_HIFI_files(sample_id, ref_title, ref_data):
    # make the ref
    id_rm_list = find_index(sample_id, ref_title)
    # print id_rm_list
    output_ref = open("_ref", "w")
    print >> output_ref, "rsID" + "\t" + list_to_line(remove_sample_from_ref(ref_title, id_rm_list)[1:])
    for pos in sorted(ref_data.keys()):
        print >> output_ref, list_to_line(remove_sample_from_ref(ref_data[pos], id_rm_list))
    output_ref.close()

    # define the person temp
    person_temp = parameter.person_dict[sample_id]

    # make the geno
    # pick_geno_list = find_parent_geno_index(sample_id, geno_title)
    # print pick_geno_list
    output_geno = open("_geno", "w")
    print >> output_geno, "rsID" + "\t" + "pos#" + "\t" + sample_id
    for pos in parameter.pos_list:
        print >> output_geno, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos]
    output_geno.close()

    # make the haplo
    # pick_haplo_list = find_parent_haplo_index(sample_id, haplo_title)
    # print pick_haplo_list
    output_haplo = open("_haplo", "w")
    print >> output_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id+"_A"
    for pos in parameter.pos_list:
        if pos in person_temp.haplotype:
            if person_temp.haplotype[pos][0] != "N" and person_temp.haplotype[pos][0] != "X":
                print >> output_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.haplotype[pos][0]
    output_haplo.close()

    # make the std_haplo
    output_std_haplo = open("_std_haplo", "w")
    print >> output_std_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id+"_A" + "\t" + sample_id+"_B"
    for pos in parameter.pos_list:
        print >> output_std_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos][0] + "\t" + person_temp.genotype_dict[pos][1]
    output_std_haplo.close()

    print sample_id, "HIFI files ready"


# def make_std_haplo_files(sample_id):
#     # define the person temp
#     person_temp = parameter.person_dict[sample_id]
#     output_std_haplo = open("_std_haplo", "w")
#     print >> output_std_haplo, "rsID" + "\t" + "pos#" + "\t" + sample_id+"_A" + "\t" + sample_id+"_B"
#     for pos in parameter.pos_list:
#         print >> output_std_haplo, parameter.rsID_dict[pos] + "\t" + str(pos) + "\t" + person_temp.genotype_dict[pos][0] + "\t" + person_temp.genotype_dict[pos][1]
#     output_std_haplo.close()


def key_compare(k1, k2):
    temp = {}
    for a in k1.keys():
        if a in k2:
            temp[a] = 0
            # print "comon key is ", a
    return temp


def cut_by_common_append():
        ref_title, ref_data = load_raw_data("refHaplos_after_merge")
        geno_title, geno_data = load_raw_data("genotype_after_merge")
        haplo_title, haplo_data = load_raw_data("haplotype_after_merge")
        com_id1 = key_compare(ref_data, geno_data)
        com_id2 = key_compare(com_id1, haplo_data)
        mini_p = min(com_id2.keys())
        max_p = max(com_id2.keys())
        # print mini_p
        # print type(mini_p)
        # print max_p
        # print type(max_p)
        geno_output = open("genotype.txt", "w")
        haplo_output = open("haplotype.txt", "w")
        ref_output = open("refHaplos.txt", "w")
        ref_list_ori = sorted(ref_data, key=int)
        geno_list_ori = sorted(geno_data, key=int)
        haplo_list_ori = sorted(haplo_data, key=int)
        ref_list = []
        geno_list = []
        haplo_list = []
        # print haplo_list
        for i in ref_list_ori:
            if int(mini_p) <= int(i) <= int(max_p):
                ref_list.append(i)
        for j in geno_list_ori:
            if int(mini_p) <= int(j) <= int(max_p):
                geno_list.append(j)
        for q in haplo_list_ori:
            if int(mini_p) <= int(q) <= int(max_p):
                haplo_list.append(q)
        # print "new list"
        # print "ref_list[0]", ref_list[0]
        # print "ref_list[-1]", ref_list[-1]
        # print "geno_list[0]", geno_list[0]
        # print "geno_list[-1]", geno_list[-1]
        # print "haplo_list[0]", haplo_list[0]
        # print "haplo_list[-1]", haplo_list[-1]
        print >> geno_output, "rsID" + "\t" + "phys_position" + "\t" + list_to_line(geno_title[2:])
        print >> haplo_output, "rsID" + "\t" + "phys_position" + "\t" + list_to_line(haplo_title[2:])
        print >> ref_output, "rsID" + "\t" + "phys_position"
        for k in sorted(ref_list):
            print >> ref_output, list_to_line(ref_data[k])
        for d in sorted(geno_list):
            print >> geno_output, list_to_line(geno_data[d])
        for c in sorted(haplo_list):
            print >> haplo_output, list_to_line(haplo_data[c])
        geno_output.close()
        haplo_output.close()
        ref_output.close()


def subprocess_execute(command, time_out=600):
    c = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    # print c.communicate()
    t = 0
    while t <= time_out and c.poll() is None:
        time.sleep(1)  # (comment 1)
        t += 1
        returncode = 1
        # print t
    if c.poll() is None:
        c.terminate()
        returncode = 0  # (comment 2)
    else:
        eturncode = c.poll()
        # print "c.poll()", c.poll()
    return returncode


def check_hapfile_run_HIFI():
    os.chdir(current_path)
    if os.path.isfile("child_hap.txt") and os.path.isfile("parent_hap.txt"):
        print "prepare for HIFI running..."
        ref_title, ref_data = load_raw_data(ref_name)
        # print "parent_hap.txt and child_hap.txt loaded..."
        sample_id_list = parameter.person_dict.keys()
        # print sample_id_list
        
        # processing samples data start
        # sample_id_list = creat_sample_id_list(parent_title)
        # print sample_id_list
        # print ref_title
        for sample_id in set(sample_id_list):
            print sample_id, "processing..."

            if os.path.exists(sample_id):  # remove the exist folder and create a new one
                shutil.rmtree(sample_id)
            os.makedirs(sample_id)

            os.chdir(sample_id)
            # print sample_id
            make_HIFI_files(sample_id, ref_title, ref_data)
            print "refMerging..."
            print subprocess.Popen("python " + refMerger + " -i _haplo -n _geno -r _ref", shell=True, stdout=subprocess.PIPE).stdout.read()

            os.rename("refHaplos.txt", "refHaplos_after_merge")
            os.rename("genotype.txt", "genotype_after_merge")
            os.rename("haplotype.txt", "haplotype_after_merge")
            cut_by_common_append()  # format required by hifi
            os.system("unix2dos refHaplos.txt")
            os.system("unix2dos genotype.txt")
            os.system("unix2dos haplotype.txt")
            hifi_starttime = time.time()
            hifi_run_code = subprocess_execute(hifi_file)
            print "hifi run time with", sample_id, round((time.time() - hifi_starttime), 6), "s"
            # print hifi_run_code
            if hifi_run_code != 0:
                print "==done=="
                print ""
            else:
                print "hifi has an issue", sample_id
                pass
            os.chdir(current_path)
        # processing samples data end
    else:
        print "Warning! child_hap.txt or parent_hap.txt file missing..."

if __name__ == '__main__':
    options = get_args()
    pedi_name = options.pedi_name
    geno_name = options.geno_file
    ref_name = options.ref_file
    global person_dict
    person_dict = {}
    
    global parameter
    parameter = parameters()
    
    start_time = time.time()
    
    genome_laser(pedi_name, geno_name)

    ###################
    # # call HIFI # #
    ###################
    print "Laser I II done: ", round(time.time() - start_time, 2), "s"
    current_path = os.getcwd()
    output_path = current_path + "/output_file"
    ref_folder = current_path + "/ref_pool"
    hifi_file = current_path + "/scripts/hifi_fu_ref_for_HIFILOCA.ref"
    refMerger = current_path + "/scripts/refMerger_v5_wli_filltheend_remove_extra.py"
    acc_check_file = current_path + "/scripts/hifiAccuCheck_v3_pos_num.py"
    # #HIFI processing
    check_hapfile_run_HIFI()
    print "elapsed_time is: ", round(time.time() - start_time, 2), "s"