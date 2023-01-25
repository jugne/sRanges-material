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


def set_topology(k_, l_, n_, j_, end_strat_ranges_, desc_on_straight_line_, parent_on_straight_line_,
                 is_branch_strat_range_, b_start_, b_end_):
    global k, l, n, j, end_strat_ranges, desc_on_straight_line, parent_on_straight_line, \
        is_branch_strat_range, b_start, b_end
    k = k_
    l = l_
    n = n_
    j = j_
    end_strat_ranges = end_strat_ranges_
    desc_on_straight_line = desc_on_straight_line_
    parent_on_straight_line = parent_on_straight_line_
    is_branch_strat_range = is_branch_strat_range_
    b_start = b_start_
    b_end = b_end_


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


def fbdr(t_or, x1=None):
    b_start_ = b_start.copy()
    b_end_ = b_end.copy()
    b_start_.insert(0, t_or)
    if x1 is not None:
        b_start_.insert(1, x1)
        b_start_.insert(2, x1)
        b_end_.insert(0, x1)
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
        if is_branch_strat_range[i]:
            second_term *= np.exp(sampling_rate * (b_start_[i] - b_end_[i]))

    return first_term * second_term * third_term * fourth_term


set_topology(k_=2, l_=1, n_=2, j_=1, end_strat_ranges_=[y1, y3],
             desc_on_straight_line_=[y2], parent_on_straight_line_=[y1], is_branch_strat_range_=[False, False, True],
             b_start_=[y1, y2], b_end_=[y1, y2, y3])

top_1 = scipy.integrate.quad(fbdr, 2., 10, args=(None,))[0]

set_topology(k_=2, l_=1, n_=2, j_=0, end_strat_ranges_=[y1, y3],
             desc_on_straight_line_=[], parent_on_straight_line_=[], is_branch_strat_range_=[False, False, False, True],
             b_start_=[y2], b_end_=[y1, y2, y3])

top_2 = scipy.integrate.dblquad(lambda x, y: fbdr(y, x), 2., 10, 2., lambda t_or: t_or)[0]

sum_top = top_1 + top_2 * 2
prob_1 = top_1 / sum_top
prob_2_3 = top_2 * 2 / sum_top
print(prob_1)
print(prob_2_3)
