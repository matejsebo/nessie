#!/usr/bin/env python

######################## 
######## SLYCE #########
######################## 
# Author: Matej Sebo   #
# Language: Python 2.7 #
# Version: 1.1         #
######################## 
# A library of useful  #
# functions to create, #
# modify, and merge    # 
# slices of lists.     #
########################


def main():
    print "Slyce tester"

    l1 = [0,1,2,3,4,5,6,7,8,9]
    
    
    len1 = len(l1)



    print "Testing forward list sanity. ----------------------"
    for begin in range(len1):
        for end in range(begin, len1):
            s = Slyce(l1, len1, begin, end, 1)
            if s.to_list() != l1[begin:end:1]:
                print "Error found in list rendering. ----------------------"
                print "Slyce params: ", l1, len1, begin, end, 1
                print s.to_list(), "!=", l1[begin:end:1]

    print "Testing reverse list sanity. ----------------------"
    for begin in range(len1):
        for end in range(-1,begin):
            s = Slyce(l1, len1, begin, end, -1)
            if s.to_list() != l1[begin:end:-1]:
                print "Error found in list rendering. ----------------------"
                print "Slyce params: ", l1, len1, begin, end, -1
                print s.to_list(), "!=", l1[begin:end:-1]

    print "Testing forward inversion (1 -> -1). ----------------------"
    for begin in range(len1):
        for end in range(begin, len1):
            s = Slyce(l1, len1, begin, end, 1)
            inv_s = s.invert()
            if inv_s.to_list() != list(reversed(l1[begin:end:1])):
                print "Error found in forward inversion."
                print "Slyce params: ", l1, len1, "(", begin, end, 1, ")"
                print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print list(reversed(l1[begin:end:1])), "!=", inv_s.to_list()

    print "Testing reverse inversion (-1 -> 1). ----------------------"
    for begin in range(len1):
        for end in range(-1,begin):
            s = Slyce(l1, len1, begin, end, -1)
            inv_s = s.invert()
            if inv_s.to_list() != list(reversed(l1[begin:end:-1])):
                print "Error found in reverse inversion."
                print "Slyce params: ", l1, len1, "(", begin, end, -1, ")"
                print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print list(reversed(l1[begin:end:-1])), "!=", inv_s.to_list()
   
    print "\nSlyce testing complete; commencing SlyceList testing."

    l2 = [11, 12, 13, 14]
    len2 = len(l2)

    print "Testing single-slyce insertion (-> | -> | ->). ----------------------"
    for pos in range(len2+1):
        for ins_list_length in range(1,5):
            ins_list = range(ins_list_length)
            final_list = l2[0:pos] + ins_list + l2[pos:]
            ins = SlyceList([Slyce(ins_list, ins_list_length)])

            s = SlyceList([Slyce(l2, len2)])
            f = SlyceList([Slyce(final_list, len(final_list))])
            test = s.insert(pos, ins)
            if test.to_list() != f.to_list():
                print "Error found in insertion."
                # print "Slyce params: ", l1, len1, "(", begin, end, -1, ")"
                # print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print f.to_list(), "!=", test.to_list()

    print "Testing single-slyce insertion (-> | <- | ->). ----------------------"
    for pos in range(len2+1):
        for ins_list_length in range(1,5):
            ins_list = range(ins_list_length)
            ins = SlyceList([Slyce(ins_list, ins_list_length)])
            final_list = l2[0:pos] + list(reversed(ins_list)) + l2[pos:]

            s = SlyceList([Slyce(l2, len2)])
            f = SlyceList([Slyce(final_list, len(final_list))])
            test = s.insert(pos, ins.invert())
            if test.to_list() != f.to_list():
                print "Error found in insertion."
                # print "Slyce params: ", l1, len1, "(", begin, end, -1, ")"
                # print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print f.to_list(), "!=", test.to_list()

    print "Testing single-slyce insertion (<- | -> | <-). ----------------------"
    for pos in range(len2+1):
        for ins_list_length in range(1,5):
            ins_list = range(ins_list_length)
            final_list = list(reversed(l2))[0:pos] + ins_list + list(reversed(l2))[pos:]
            ins = SlyceList([Slyce(ins_list, ins_list_length)])

            s = SlyceList([Slyce(l2, len2)]).invert()
            f = SlyceList([Slyce(final_list, len(final_list))])
            test = s.insert(pos, ins)
            if test.to_list() != f.to_list():
                print "Error found in insertion."
                # print "Slyce params: ", l1, len1, "(", begin, end, -1, ")"
                # print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print f.to_list(), "!=", test.to_list()

    print "Testing single-slyce insertion (<- | <- | <-). ----------------------"
    for pos in range(len2+1):
        for ins_list_length in range(1,5):
            ins_list = range(ins_list_length)
            ins = SlyceList([Slyce(ins_list, ins_list_length)])
            final_list = list(reversed(l2))[0:pos] + list(reversed(ins_list)) + list(reversed(l2))[pos:]

            s = SlyceList([Slyce(l2, len2)]).invert()
            f = SlyceList([Slyce(final_list, len(final_list))])
            test = s.insert(pos, ins.invert())
            if test.to_list() != f.to_list():
                print "Error found in insertion."
                # print "Slyce params: ", l1, len1, "(", begin, end, -1, ")"
                # print "Invert params: ", inv_s.l, inv_s.length, "(", inv_s.i1, inv_s.i2, inv_s.i3, ")"
                print f.to_list(), "!=", test.to_list()

    l3 = [30,31,32,33,34,35,36,37,38,39]
    len3 = len(l3)

    print "Testing single-slyce forward extraction. ----------------------"
    for ext_list_length in range(len3):
        for pos in range(len3-ext_list_length+1):
            ext_list = l3[pos:pos+ext_list_length]
            s = SlyceList([Slyce(l3, len3)])
            f = SlyceList([Slyce(ext_list, len(ext_list))])
            test = s.extract(pos, ext_list_length)
            if test.to_list() != f.to_list():
                print "Error found in extraction."
                print f.to_list(), "!=", test.to_list()

    print "Testing single-slyce reverse extraction. ----------------------"
    for ext_list_length in range(len3):
        for pos in range(len3-ext_list_length):
            ext_list = list(reversed(l3))[pos:pos+ext_list_length]
            s = SlyceList([Slyce(l3, len3)])
            s = s.invert()
            f = SlyceList([Slyce(ext_list, len(ext_list))])
            test = s.extract(pos, ext_list_length)
            if test.to_list() != f.to_list():
                print "Error found in extraction."
                print f.to_list(), "!=", test.to_list()


    print "Testing single-slyce forward excision. ----------------------"
    #print l3
    for ext_list_length in range(len3):
        for pos in range(len3-ext_list_length+1):
            #print "Excising at", pos, "with length", ext_list_length
            ext_list = l3[pos:pos+ext_list_length]
            np_list = l3[:pos] + l3[pos+ext_list_length:]
            s = SlyceList([Slyce(l3, len3)])
            f = SlyceList([Slyce(ext_list, len(ext_list))])
            f_np = SlyceList([Slyce(np_list, len(np_list))])
            np, extr = s.excise(pos, ext_list_length)
            if extr.to_list() != f.to_list():
                print "Error found in excision extract."
                print f.to_list(), "!=", extr.to_list()
            if np.to_list() != f_np.to_list():
                print "Error found in excision new parent."
                print f_np.to_list(), "!=", np.to_list()

    print "Testing single-slyce reverse excision. ----------------------"
    for ext_list_length in range(len3):
        for pos in range(len3-ext_list_length):
            ext_list = list(reversed(l3))[pos:pos+ext_list_length]
            np_list = list(reversed(l3))[:pos] + list(reversed(l3))[pos+ext_list_length:]
            s = SlyceList([Slyce(l3, len3)])
            s = s.invert()
            f = SlyceList([Slyce(ext_list, len(ext_list))])
            f_np = SlyceList([Slyce(np_list, len(np_list))])
            np, extr = s.excise(pos, ext_list_length)
            if extr.to_list() != f.to_list():
                print "Error found in excision extract."
                print f.to_list(), "!=", extr.to_list()
            if np.to_list() != f_np.to_list():
                print "Error found in excision new parent."
                print f_np.to_list(), "!=", np.to_list()

    l4 = [40,41,42,43]
    len4 = len(l4)
    l5 = [50,51,52,53]
    len5 = len(l5)
    l6 = [60,61,62,63]
    len6 = len(l6)
    l_tot = l4+l5+l6
    len_tot = len(l_tot)
    sl_tot = SlyceList([Slyce(l4, len4), Slyce(l5, len5), Slyce(l6, len6)])
    
    print "Testing multi-slyce forward excision (-|> -> -|>). ----------------------"
    #print l3
    for pos in range(0, 3):
        for len_extr in range(4, 12):
            # print "Excising at", pos, "with length", len_extr
            ext_list = l_tot[pos:pos+len_extr]
            np_list = l_tot[:pos] + l_tot[pos+len_extr:]

            np, extr = sl_tot.excise(pos, len_extr)
            if extr.to_list() != ext_list:
                print "Error found in excision extract."
                print ext_list, "!=", extr.to_list()
            if np.to_list() != np_list:
                print "Error found in excision new parent."
                print np_list, "!=", np.to_list()

    l4 = [40,41,42,43]
    len4 = len(l4)
    l5 = [50,51,52,53]
    len5 = len(l5)
    l6 = [60,61,62,63]
    len6 = len(l6)
    l_tot = list(reversed(l4))+list(reversed(l5))+list(reversed(l6))
    len_tot = len(l_tot)
    sl_tot = SlyceList([Slyce(l4, len4).invert(), Slyce(l5, len5).invert(), Slyce(l6, len6).invert()])

    print "Testing multi-slyce reverse excision (<|- <- <|-). ----------------------"
    #print l3
    for pos in range(0, 3):
        for len_extr in range(4, 12):
            # print "Excising at", pos, "with length", len_extr
            ext_list = l_tot[pos:pos+len_extr]
            np_list = l_tot[:pos] + l_tot[pos+len_extr:]

            np, extr = sl_tot.excise(pos, len_extr)
            if extr.to_list() != ext_list:
                print "Error found in excision extract."
                print ext_list, "!=", extr.to_list()
            if np.to_list() != np_list:
                print "Error found in excision new parent."
                print np_list, "!=", np.to_list()

    l3 = [30,31,32,33]
    len3 = len(l3)
    l4 = [40,41,42,43]
    len4 = len(l4)
    l5 = [50,51,52,53]
    len5 = len(l5)
    l6 = [60,61,62,63]
    len6 = len(l6)
    l7 = [70,71,72,73]
    len7 = len(l7)
    l_tot = l3+l4+l5+l6+l7
    len_tot = len(l_tot)
    sl_tot = SlyceList([Slyce(l3, len3), Slyce(l4, len4), Slyce(l5, len5), Slyce(l6, len6), Slyce(l7, len7)])
    
    print "Testing multi-slyce forward excision (-> -|> -> -|> ->). ----------------------"
    #print l3
    for pos in range(4,7):
        for len_extr in range(4, 8):
            # print "Excising at", pos, "with length", len_extr
            ext_list = l_tot[pos:pos+len_extr]
            np_list = l_tot[:pos] + l_tot[pos+len_extr:]

            np, extr = sl_tot.excise(pos, len_extr)
            if extr.to_list() != ext_list:
                print "Error found in excision extract."
                print ext_list, "!=", extr.to_list()
            if np.to_list() != np_list:
                print "Error found in excision new parent."
                print np_list, "!=", np.to_list()

    l3 = [30,31,32,33]
    len3 = len(l3)
    l4 = [40,41,42,43]
    len4 = len(l4)
    l5 = [50,51,52,53]
    len5 = len(l5)
    l6 = [60,61,62,63]
    len6 = len(l6)
    l7 = [70,71,72,73]
    len7 = len(l7)
    l_tot = list(reversed(l3))+list(reversed(l4))+list(reversed(l5))+list(reversed(l6))+list(reversed(l7))
    len_tot = len(l_tot)
    sl_tot = SlyceList([Slyce(l3, len3).invert(), Slyce(l4, len4).invert(), Slyce(l5, len5).invert(), Slyce(l6, len6).invert(), Slyce(l7, len7).invert()])

    print "Testing multi-slyce reverse excision (<- <|- <- <|- <-). ----------------------"
    #print l3
    for pos in range(4, 7):
        for len_extr in range(4, 8):
            # print "Excising at", pos, "with length", len_extr
            ext_list = l_tot[pos:pos+len_extr]
            np_list = l_tot[:pos] + l_tot[pos+len_extr:]

            np, extr = sl_tot.excise(pos, len_extr)
            if extr.to_list() != ext_list:
                print "Excising at", pos, "with length", len_extr
                print "Error found in excision extract."
                print ext_list, "!=", extr.to_list()
            if np.to_list() != np_list:
                print "Excising at", pos, "with length", len_extr
                print "Error found in excision new parent."
                print np_list, "!=", np.to_list()

# These slyce function assume:
#   - Non-negative slyce index i1
#   - Positive, 0, -1, or None slyce index i2
#   - slyce index i3 that is 1 or -1
#   - A scaffold id "i" corresponding to the metadata of a particular scaffold
#       - In the case of the dump function (for debugging), i is the parent list

def len_scaf(scaf_id):
    return len(scaf_id)

class Slyce(object):
    def __init__(self, l, length, i1=None, i2=None, i3=1, p=0):
        self.l = l # parent list identifier

        if i3 == None or i3 > 0:
            i3 = 1
        else:
            i3 = -1

        if i1 == None or i2 == None:
            self.i1, self.i2, self.i3 = \
                slice(i1, i2, i3).indices(length)
        else:
            self.i1, self.i2, self.i3 = i1, i2, i3
            
        self.length = length # length of the parent list
        self.p = p # partition

    @staticmethod
    def force_init(l, length, i1, i2, i3, p=0):
        s = Slyce(l, length)
        s.i1 = i1
        s.i2 = i2
        s.i3 = i3
        s.p = p
        return s

    @staticmethod
    def blank_slice(l, length, p=0):
        return Slyce(l, length, 0, 0, 1, p)

    def invert(self):
        if self.i2 == self.i1: return Slyce.blank_slice(self.l, self.length, self.p)
        new_i1, new_i2, new_i3 = None, None, None
        if self.i3 == 1:
            new_i1 = None if self.i2 == None or self.i2 < 0 else self.i2-1
            new_i2 = None if self.i1 == None or self.i1 < 0 else self.i1-1
            new_i3 = -1
        else: # if i3 == -1:
            new_i1 = None if self.i2 == None else self.i2+1
            new_i2 = None if self.i1 == None else self.i1+1
            new_i3 = 1

        if self.i2 == -1 and self.i1 != -1 and new_i3 == 1: return Slyce.blank_slice(self.l, self.length, self.p)
        if new_i2 == -1 and new_i1 != -1 and self.i3 == 1: return Slyce.force_init(self.l, self.length, new_i1, None, new_i3, self.p)
        return Slyce(self.l, self.length, new_i1, new_i2, new_i3, self.p)

    def all_after_pos(self, pos): # INCLUSIVE W. RESPECT TO POS
        # print pos, self.i1, self.i2
        if self.i3 == 1:
            new_i1 = self.i1 + pos
            new_i2 = self.i2
        else: 
            new_i1 = self.i1 - pos
            new_i2 = self.i2
        if new_i2 == None or new_i2 <= -1:
            return Slyce.force_init(self.l, self.length, new_i1, None, self.i3, self.p)
        return Slyce(self.l, self.length, new_i1, new_i2, self.i3, self.p)

    def all_before_pos(self, pos): # NON-INCLUSIVE W. RESPECT TO POS
        # print pos, self.i1, self.i2
        if self.i3 == 1:
            new_i1 = self.i1
            new_i2 = self.i1+pos
        else: 
            new_i1 = self.i1
            new_i2 = self.i1 - pos 
        if new_i2 == None or new_i2 <= -1:
            return Slyce.force_init(self.l, self.length, new_i1, None, self.i3, self.p)
        return Slyce(self.l, self.length, new_i1, new_i2, self.i3, self.p)

    # This function assumes pos1 < pos2
    def all_between_pos(self, pos1, pos2): # INCLUSIVE W. RESPECT TO POS1, NON-INCLUSIVE W. RESPECT TO POS2
        # print pos, self.i1, self.i2
        if self.i3 == 1:
            new_i1 = self.i1+pos1
            new_i2 = self.i1+pos2
        else: # if i3 == -1:# following code is inclusive w respect to pos
            new_i1 = self.i1-pos1
            new_i2 = self.i1-pos2
        if new_i2 == None or new_i2 <= -1:
            return Slyce.force_init(self.l, self.length, new_i1, None, self.i3, self.p)
        return Slyce(self.l, self.length, new_i1, new_i2, self.i3, self.p)

    def len(self):
        return abs(self.i1-(-1 if self.i2 == None else self.i2))

    def to_slice(self):
        i1, i2, i3 = self.i1, self.i2, self.i3
        slyce_i = (i1, i2, i3)
        return slice(*slyce_i)

    def post_continuous(self, other_slyce):
        if self.l == other_slyce.l and self.i3 == other_slyce.i3 and \
            self.i2 == other_slyce.i1:
                return True
        return False
    def continuous(self, other_slyce):
        if self.l == other_slyce.l and self.i3 == other_slyce.i3 and \
            (self.i2 == other_slyce.i1 or self.i1 == other_slyce.i2):
                return True
        return False

    # increase length of slyce by one (via postfix)
    def add_one(self):
        if self.i3 == 1:
            self.i2 += 1
        else: # if self.i3 == -1:
            self.i1 -= 1

    # decrease length of slyce by one (via postdel)
    def zap_one(self):
        if self.i3 == 1:
            self.i2 -= 1
        else: # if self.i3 == -1:
            self.i1 += 1

    # The following function (used for debugging only!) assumes l in (l,len,i1,i2,i3) is a valid parent list
    def to_list(self):
        return self.l[self.to_slice()]

    def __str__(self):
        return "Slyce" + str((self.l, self.length, (self.i1, self.i2, self.i3), self.p))

class SlyceList(object):
    def __init__(self, l=None):
        if not l:
            self.sl = []
        else:
            self.sl = l
        
    def invert(self):
        return SlyceList(list(reversed([s.invert() for s in self.sl])))

    def len(self):
        return sum([s.len() for s in self.sl])

    def cat(self, other):
        self.sl += other.sl

    def slyce_index_at_abs_pos(self, pos):
        total_bases = 0
        slyce_index = 0
        while total_bases < pos and slyce_index < len(self.sl):
            total_bases += self.sl[slyce_index].len()
            slyce_index += 1
        slyce_index -= 1
        if slyce_index == -1:
            slyce_index = 0
            total_bases = self.sl[slyce_index].len()
        return slyce_index, total_bases

    # insert slyce list ("ins") into self at absolute position "position"
    def insert(self, position, ins):
        # print "INSERT", self.to_list(), position, ins.to_list()
        slyce_index, total_bases = self.slyce_index_at_abs_pos(position)
        last_slyce_len = self.sl[slyce_index].len()

        rel_loc = position - total_bases + last_slyce_len
        new_slyce_list_here = \
            [self.sl[slyce_index].all_before_pos(rel_loc)] + \
            ins.sl + \
            [self.sl[slyce_index].all_after_pos(rel_loc)]

        return SlyceList(self.sl[0:slyce_index] + new_slyce_list_here + self.sl[max(1,slyce_index+1):]).pare()

    # does not modify self, returns extracted SlyceList
    def extract(self, position, ex_len):
        slyce_index_start, total_bases_start = self.slyce_index_at_abs_pos(position)
        slyce_index_end, total_bases_end = self.slyce_index_at_abs_pos(position + ex_len)

        if slyce_index_start == slyce_index_end and slyce_index_start:
            rel_pos = position - total_bases_start + self.sl[slyce_index_start].len()
            extract = SlyceList([self.sl[slyce_index_start].all_between_pos(rel_pos, rel_pos+ex_len)])
        else:
            #print "b"
            start_slyce_len = self.sl[slyce_index_start].len()
            end_slyce_len = self.sl[slyce_index_end].len()

            rel_pos_start = position - total_bases_start + start_slyce_len
            rel_pos_end = position + ex_len - total_bases_end + end_slyce_len

            if slyce_index_start == slyce_index_end:
                extract = SlyceList([self.sl[slyce_index_start].all_between_pos(rel_pos_start, rel_pos_end)])
            else: 
                extract = SlyceList([self.sl[slyce_index_start].all_after_pos(rel_pos_start)] + \
                    self.sl[slyce_index_start+1:slyce_index_end] + \
                    [self.sl[slyce_index_end].all_before_pos(rel_pos_end)])
        return extract.pare()

    # does not modify self, returns original SlyceList without the extract AND the extracted SlyceList
    def excise(self, position, ex_len):
        #print "EXCISE", self.to_list(), position, ex_len
        slyce_index_start, total_bases_start = self.slyce_index_at_abs_pos(position)
        slyce_index_end, total_bases_end = self.slyce_index_at_abs_pos(position + ex_len)

        if slyce_index_start == slyce_index_end and slyce_index_start:
            rel_pos = position - total_bases_start + self.sl[slyce_index_start].len()
            extract = SlyceList([self.sl[slyce_index_start].all_between_pos(rel_pos, rel_pos+ex_len)])
            new_parent = SlyceList(self.sl[0:slyce_index_start] + \
                [self.sl[slyce_index_start].all_before_pos(rel_pos)] + \
                [self.sl[slyce_index_start].all_after_pos(rel_pos+ex_len)] + \
                self.sl[slyce_index_end+1:])
        else:
            start_slyce_len = self.sl[slyce_index_start].len()
            end_slyce_len = self.sl[slyce_index_end].len()

            rel_pos_start = position - total_bases_start + start_slyce_len
            rel_pos_end = position + ex_len - total_bases_end + end_slyce_len

            if slyce_index_start == slyce_index_end:
                extract = SlyceList([self.sl[slyce_index_start].all_between_pos(rel_pos_start, rel_pos_end)])
                new_parent = SlyceList(self.sl[0:slyce_index_start] + \
                    [self.sl[slyce_index_start].all_before_pos(rel_pos_start)] + \
                    [self.sl[slyce_index_end].all_after_pos(rel_pos_end)] + \
                    self.sl[slyce_index_end+1:])
            else: 
                extract = SlyceList([self.sl[slyce_index_start].all_after_pos(rel_pos_start)] + \
                    self.sl[slyce_index_start+1:slyce_index_end] + \
                    [self.sl[slyce_index_end].all_before_pos(rel_pos_end)])
                new_parent = SlyceList(self.sl[0:slyce_index_start] + \
                    [self.sl[slyce_index_start].all_before_pos(rel_pos_start)] + \
                    [self.sl[slyce_index_end].all_after_pos(rel_pos_end)] + \
                    self.sl[slyce_index_end+1:])
        return new_parent.pare(), extract.pare()

    # remove spurious (empty/malformed) slyces from a slyce list
    def pare(self):
        new_s = []
        for s in self.sl:
            if s.i1 != s.i2 and (s.i1 == None or s.i1 > -1) and (s.i2 == None or s.i2 >= -1):
                new_s += [s]
        return SlyceList(new_s)

    # If back-to-back slyces are contiguous and the SlyceList can be shortened, do so (and modify self).
    def join_ends(self):
        if not self.sl:
            return
        new_sl = [self.sl[0]]
        top_index = 0
        for added_slyce in range(1,len(self.sl)):
            if new_sl[top_index].post_continuous(added_slyce):
                new_sl[top_index].i2 = added_slyce
            else:
                new_sl += [added_slyce]
                top_index += 1
        self.sl = new_sl

    # The length of all slyces inside this SlyceList that belong to partition p.
    def len_partition(self, p):
        total = 0
        for s in self.sl:
            if s.p == p:
                total += s.len()
        return total

    # The following functions (used ONLY for debugging) assume l in (l,i1,i2,i3) is a valid parent list
    def to_list(self):
        all_slyces = []
        for s in self.sl:
            all_slyces += s.to_list()
        return all_slyces

    # print evaluated contents of a slyce list
    def dump(self):
        print self.to_list()

    # print the unevaluated representation of a slyce list
    def __str__(self):
        if not self.sl:
            return "SlyceList [<BLANK>]"
        r = "SlyceList[\n"
        for s in self.sl:
            r += "   " + str(s) + "\n"
        return r + "]"

if __name__ == '__main__':
    main()
