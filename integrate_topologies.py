import csv
import numpy as np
import scipy.integrate

birth_rate = 2.0
death_rate = 1.0
sampling_rate = 0.5
removal_rate = 0.9
sampling_present_prob = 0.5

# total number of sampled fossils
k = 2
# Number of sampled stratigraphic ranges with an extant species sample
l = 1
# Number of sampled stratigraphic ranges, i.e. , number of sampled species
# (some stratigraphic ranges may only be represented by a single sample
# in the past or present)
n = 2
# Number of sampled-ancestor-stratigraphic ranges (the whole range is SA to other range)
j = 1

y1 = 2.0
y2 = 1.0
y3 = 0.0
end_strat_ranges = []
# Time of first observed fossil corresponding to the species represented by
# stratigraphic range i , i.e. , start time of stratigraphic range i
desc_on_straight_line = []
# Time of last observed fossil corresponding to the species represented by
# stratigraphic range i , i.e. , end time of stratigraphic range i
parent_on_straight_line = []
is_branch_strat_range = []

# Start and end of each branch
b_start = []
b_end = []

# Start and end of each branch
b_start_range = []
b_end_range = []


def set_topology(k_, l_, n_, j_, end_strat_ranges_, desc_on_straight_line_, parent_on_straight_line_,
                 is_branch_strat_range_, b_start_range_=[], b_end_range_=[]):
    global k, l, n, j, end_strat_ranges, desc_on_straight_line, parent_on_straight_line, \
        is_branch_strat_range, b_start_range, b_end_range
    k = k_
    l = l_
    n = n_
    j = j_
    end_strat_ranges = end_strat_ranges_
    desc_on_straight_line = desc_on_straight_line_
    parent_on_straight_line = parent_on_straight_line_
    is_branch_strat_range = is_branch_strat_range_
    b_start_range = b_start_range_
    b_end_range = b_end_range_


rate_diff = birth_rate - death_rate - sampling_rate
rate_sum = birth_rate + death_rate + sampling_rate

c1 = np.sqrt(rate_diff * rate_diff + 4 * birth_rate * sampling_rate)
c2 = -((rate_diff - 2 * birth_rate * sampling_present_prob) / c1)


def q(x):
    return (4 * np.exp(-c1 * x)) / \
           ((np.exp(-c1 * x) * (1 - c2) + (1 + c2)) *
            (np.exp(-c1 * x) * (1 - c2) + (1 + c2)))


def qasym(x):
    return np.sqrt(np.exp(-x * rate_sum) * q(x))


def qasym_hat(start, end, strat_range):
    if strat_range:
        return qasym(start) / qasym(end)
    else:
        return q(start) / q(end)


def p0s(x):
    return removal_rate + (1 - removal_rate) * p(x)


def p(x):
    return 1 + (-rate_diff + c1 * (
            (np.exp(-c1 * x) * (1 - c2) - (1 + c2)) /
            (np.exp(-c1 * x) * (1 - c2) + (1 + c2)))) / (2.0 * birth_rate)


def fbdr_x1(t_or, x1, topology):
    return fbdr(t_or, x1, None, topology)


def fbdr_x2(t_or, x2, topology):
    return fbdr(t_or, None, x2, topology)


def fbdr(t_or, x1=None, x2=None, topology=None):
    if topology is not None:
        if topology == 318:
            b_start_ = [t_or, y1, y2]
            b_end_ = [y1, y2, y3]
        elif topology == 317:
            b_start_ = [t_or, x1, x1, y2]
            b_end_ = [x1, y1, y2, y3]
        ######
        elif topology == 321:
            b_start_ = [t_or, x1, x1, x2, x2]
            b_end_ = [x1, y1, x2, y2, y3]
        elif topology == 322:
            b_start_ = [t_or, x1, x1, x2, x2]
            b_end_ = [x1, x2, y3, y1, y2]
        elif topology == 323:
            b_start_ = [t_or, x1, x1, x2, x2]
            b_end_ = [x1, x2, y2, y1, y3]
        elif topology == 324:
            b_start_ = [t_or, y1, x2, x2]
            b_end_ = [y1, x2, y2, y3]
        elif topology == 325:
            b_start_ = [t_or, x1, x1, y1]
            b_end_ = [x1, y1, y3, y2]
        elif topology == 326:
            b_start_ = [t_or, x1, x1, y1]
            b_end_ = [x1, y1, y2, y3]
        elif topology == 327:
            b_start_ = [t_or, x1, x1, y2]
            b_end_ = [x1, y1, y2, y3]
        elif topology == 328:
            b_start_ = [t_or, y1, y2]
            b_end_ = [y1, y2, y3]
        ########
        elif topology == 337:
            b_start_ = [t_or, y1, y2]
            b_end_ = [y1, y2, y3]
        elif topology == 333:
            b_start_ = [t_or, y1, x2, x2]
            b_end_ = [y1, x2, y2, y3]
        elif topology == 334:
            b_start_ = [t_or, x1, x1, y1]
            b_end_ = [x1, y1, y3, y2]
    # first_term = 1/(1 - p(t_or))
    first_term = np.power(sampling_rate, k) * np.power(sampling_present_prob, l) * \
                 np.power(birth_rate, (n - j - 1)) / (1 - p(t_or))

    # y_i - Time of last observed fossil corresponding to the
    # species represented by stratigraphic range i , i.e. ,
    # end time of stratigraphic range i
    third_term = 1
    for i in range(0, n - j - l):
        third_term *= p0s(end_strat_ranges[i])

    # Let i ∈ I if stratigraphic range i and its most recent ancestral stratigraphic
    # range, a ( i ), lie on a straight line in the graphical representation of
    # the sampled tree.
    fourth_term = 1
    for i in range(len(desc_on_straight_line)):
        o_i = desc_on_straight_line[i]
        y_ai = parent_on_straight_line[i]
        fourth_term *= (1 - (q(o_i) / qasym(o_i)) * (qasym(y_ai) / q(y_ai)))

    second_term = 1
    for i in range(len(b_start_)):
        second_term *= qasym_hat(b_start_[i], b_end_[i], is_branch_strat_range[i])
    for i in range(len(b_start_range)):
        second_term *= np.exp(sampling_rate * (b_start_range[i] - b_end_range[i]))

    return first_term * second_term * third_term * fourth_term

print("Topologies for 3 samples (ignoring ranges and orientation):\n "
      "1-((3,2),1); 2-((3,1),2); 3-(3,(2,1)); 4-((3,2)1); 5-(3,(2)1); 6-(2,(3)1); 7-(1,(3)2); 8-(((3)2)1).")

set_topology(k_=2, l_=1, n_=2, j_=1, end_strat_ranges_=[y1],
             desc_on_straight_line_=[y2], parent_on_straight_line_=[y1], is_branch_strat_range_=[False, False, True],
             b_start_range_=[y2], b_end_range_=[y3])

top_318 = scipy.integrate.quad(fbdr, y1, 10, args=(None, None, 318))[0]

set_topology(k_=2, l_=1, n_=2, j_=0, end_strat_ranges_=[y1],
             desc_on_straight_line_=[], parent_on_straight_line_=[], is_branch_strat_range_=[False, False, False, True],
             b_start_range_=[y2], b_end_range_=[y3])

top_317 = scipy.integrate.dblquad(fbdr_x1, y1, 10, lambda x1: x1, lambda x1: 10,
                                args=(317,))[0]

sum_top = top_318 + top_317 * 2
prob_1 = top_318 / sum_top
prob_2_3 = top_317 * 2 / sum_top
print("Probabilities for three samples-one range between the lower two samples y2,y3:")
print("Probability of topology 7 in two orientations: {}".format(prob_2_3))
print("Probability of topology 8: {}".format(prob_1))


##### Three samples-three ranges topologies


### topology 1, 4 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=0, end_strat_ranges_=[y1, y2],
             desc_on_straight_line_=[], parent_on_straight_line_=[],
             is_branch_strat_range_=[False, False, False, False, False])

# logP = np.log(fbdr(4., 2.5044396886075146, 1.6667133337533684, 1))
top_3_1 = scipy.integrate.tplquad(fbdr, y2, 10,
                                  lambda x2: max(x2, y1), lambda x2: 10,
                                  lambda x2, x1: x1, lambda x2, x1: 10,
                                  args=(321,))[0] * 4


### topology 2, 4 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=0, end_strat_ranges_=[y1, y2],
             desc_on_straight_line_=[], parent_on_straight_line_=[],
             is_branch_strat_range_=[False, False, False, False, False])
top_3_2 = scipy.integrate.tplquad(fbdr, y1, 10,
                                  lambda x2: x2, lambda x2: 10,
                                  lambda x2, x1: x1, lambda x2, x1: 10,
                                  args=(322,))[0] * 4

### topology 3, 4 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=0, end_strat_ranges_=[y1, y2],
             desc_on_straight_line_=[], parent_on_straight_line_=[],
             is_branch_strat_range_=[False, False, False, False, False])

# logP = np.log(fbdr(4.338697113659037, 2.5044396886075146, 2.1901560384944956, 5))
top_3_3 = scipy.integrate.tplquad(fbdr, y1, 10,
                                  lambda x2: x2, lambda x2: 10,
                                  lambda x2, x1: x1, lambda x2, x1: 10,
                                  args=(323,))[0] * 4

### topology 4, 2 UNEQUAL orientations
set_topology(k_=2, l_=1, n_=3, j_=1, end_strat_ranges_=[y2],
             desc_on_straight_line_=[y2], parent_on_straight_line_=[y1],
             is_branch_strat_range_=[False, False, False, False])
# logP = np.log(fbdr_x2(3.682852107065018, 1.5807170686314225, 2))
top_3_41 = scipy.integrate.dblquad(fbdr_x2, y2, y1, y1, 10,
                                   args=(324,))[0]

set_topology(k_=2, l_=1, n_=3, j_=1, end_strat_ranges_=[y2],
             desc_on_straight_line_=[y3], parent_on_straight_line_=[y1],
             is_branch_strat_range_=[False, False, False, False])
# logP = np.log(fbdr_x2(3.682852107065018, 1.5807170686314225, 2))
top_3_42 = scipy.integrate.dblquad(fbdr_x2, y2, y1, y1, 10,
                                   args=(324,))[0]

### topology 5, 2 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=1, end_strat_ranges_=[y2],
             desc_on_straight_line_=[y2], parent_on_straight_line_=[y1],
             is_branch_strat_range_=[False, False, False, False])
top_3_5 = scipy.integrate.dblquad(fbdr_x1, y1, 10, lambda x1: x1, lambda x1: 10,
                                  args=(325,))[0] * 2

### topology 6, 2 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=1, end_strat_ranges_=[y2],
             desc_on_straight_line_=[y3], parent_on_straight_line_=[y1],
             is_branch_strat_range_=[False, False, False, False])
top_3_6 = scipy.integrate.dblquad(fbdr_x1, y1, 10, lambda x1: x1, lambda x1: 10,
                                  args=(326,))[0] * 2

### topology 7, 2 equal orientations
set_topology(k_=2, l_=1, n_=3, j_=1, end_strat_ranges_=[y2],
             desc_on_straight_line_=[y3], parent_on_straight_line_=[y2],
             is_branch_strat_range_=[False, False, False, False])
top_3_7 = scipy.integrate.dblquad(fbdr_x1, y1, 10, lambda x1: x1, lambda x1: 10,
                                  args=(327,))[0] * 2

### topology 8, only one orientation
set_topology(k_=2, l_=1, n_=3, j_=2, end_strat_ranges_=[],
             desc_on_straight_line_=[y2, y3], parent_on_straight_line_=[y1, y2],
             is_branch_strat_range_=[False, False, False])
top_3_8 = scipy.integrate.quad(fbdr, 2., 10, args=(None, None, 328))[0]

sum_ = top_3_1 + top_3_2 + top_3_3 + top_3_41 + top_3_42 + top_3_5 + top_3_6 + top_3_7 + top_3_8
print("Probabilities for three samples-three ranges:")
print("Probability of topology 1 in 4 (equal) orientations: {}".format(top_3_1 / sum_))
print("Probability of topology 2 in 4 (equal) orientations: {}".format(top_3_2 / sum_))
print("Probability of topology 3 in 4 (equal) orientations: {}".format(top_3_3 / sum_))
print("Probability of topology 4 in 2 (unequal) orientations: {0}={1}+{2}".format((top_3_41+top_3_42) / sum_,
                                                                                  top_3_41/sum_, top_3_42/sum_))
print("Probability of topology 5 in 2 (equal) orientations: {}".format(top_3_5 / sum_))
print("Probability of topology 6 in 2 (equal) orientations: {}".format(top_3_6 / sum_))
print("Probability of topology 7 in 2 (equal) orientations: {}".format(top_3_7 / sum_))
print("Probability of topology 8: {}".format(top_3_8 / sum_))




### topology 8, only one orientation
set_topology(k_=2, l_=1, n_=2, j_=1, end_strat_ranges_=[],
             desc_on_straight_line_=[y3], parent_on_straight_line_=[y2],
             is_branch_strat_range_=[False, True, False],
             b_start_range_=[y1], b_end_range_=[y2])
top_32_8 = scipy.integrate.quad(fbdr, 2., 10, args=(None, None, 328))[0]

### topology 4, only one orientation
set_topology(k_=2, l_=1, n_=2, j_=0, end_strat_ranges_=[y2],
             desc_on_straight_line_=[], parent_on_straight_line_=[],
             is_branch_strat_range_=[False, True, True, False],
            b_start_range_=[y1], b_end_range_=[y2])
top_32_4 = scipy.integrate.dblquad(fbdr_x2, y2, y1, y1, 10, args=(324,))[0]

### topology 5, only one orientation
set_topology(k_=2, l_=1, n_=2, j_=0, end_strat_ranges_=[y2],
             desc_on_straight_line_=[], parent_on_straight_line_=[],
             is_branch_strat_range_=[False, False, False, True],
             b_start_range_=[y1], b_end_range_=[y2])
top_32_5 = scipy.integrate.dblquad(fbdr_x1, y1, 10, lambda x1: x1, lambda x1: 10,
                                args=(325,))[0]*2

sum_ = top_32_8 + top_32_4 + top_32_5
print("Probabilities for three samples-one range between the top two samples y1,y2:")
print("Probability of topology 4: {}".format(top_32_4 / sum_))
print("Probability of topology 5 in 2 (equal) orientations: {}".format(top_32_5 / sum_))
print("Probability of topology 8: {}".format(top_32_8 / sum_))
