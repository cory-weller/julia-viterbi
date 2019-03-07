# Julia Implementation of Viterbi Algorithm - jho5ze 1st draft 3-6-19

#Citations:
#    https://en.wikipedia.org/wiki/Viterbi_algorithm
#    https://www.robots.ox.ac.uk/~vgg/rg/papers/hmm.pdf
#    https://scls.gitbooks.io/ljthw/content/_chapters/04-ex1.html
#    https://en.wikibooks.org/wiki/Introducing_Julia/Working_with_text_files


#USE ARG_MAX INSTEAD OF ITERATIONS FOR FASTER IMPLEMENTATION

#ARGS holds the command line arguments - can be used to input a file to be read in

#NEED TO INITIALIZE VARIABLES WITH PROPER TYPES AND SIZES

#O = observation space - sequence of nucleotides (or SNP variants) at a certain location - may be one or two SNPs based on if get read from both chromosomes or just one
# S = state space - haplotypes
# pi = initial probabilities of the state spaces - probability that the initial state starts at the given state (1D array with numbered haplotypes being the index and prob being at that index) - ACTUALLY just gonna say theyre 1 since all haplotypes equally likely AT THE MOMENT
# A = transition matrix - probability that the state will transition between one state to the next - Only possible to stay in the same two haplotype state or change one of two since recombination of both at the same site (or between the SNPs being assessed is unlikely) - will just be a placeholder 0.9 for staying in the same state and 0.1/#of other possibilities for the rest OR just 0.1 for general crossover (NEED TO DECIDE HERE)- future implementations can also factor in the %recombination between the SNPs - generated on the fly based on these values (not an actual matrix - implement as a function or statically)
# B = emission matrix - probability of observing the observation from a certain state - i.e. getting a certain haplotype you get based on the SNP you read - TENTATIVE: 0.45 for one chromosome having proper allele, 0.9 for both, 0.1 for neither - Also generated on the fly - implement as a function
#T1 = stores the probability of the most likely path so far with the last element being that designated by its row (which is a specific haplotype). Using the convention that the T1 and T2 tables are being filled in with the first label column being AA corresponding to 1,AB corresponding to 2,AC corresponding to 3,AD,...,BA,BB,BC,... etc.
#T2 = stores the last element of the most likely path so far to be tarversed using the max of the column from above to find the ultimate solution. Using the convention that the T1 and T2 tables are being filled in with the first label column being AA corresponding to 1,AB corresponding to 2,AC corresponding to 3,AD,...,BA,BB,BC,... etc.
#LINE = current line of the file being read - observation at index obs, haplotype data starts at index obs + 1, and there are num_haplotypes haplotypes


println(string(ARGS[1]))

observation_file = open(string(ARGS[1]))
haplotypes_file = open(string(ARGS[2]))

observations = readlines(observation_file)
haplotypes = readlines(haplotypes_file)

obs_index = 2
global hap_index = 5

num_SNPs = length(haplotypes) #(number of data lines in the file)

header_names = split(haplotypes[1])

line_counter = 2
global HAPS_LINE = split(haplotypes[2])
OBS_ALLELE = (split(observations[2])[1]
global num_haplotypes = length(HAPS_LINE)

T1 = Array{Float64}(undef, num_haplotypes^2, num_SNPs - 1)
T2 = Array{Float64}(undef, num_haplotypes^2, num_SNPs - 1)

#Initialize the first value in T1 (the initial probabilities (all just 1 for now) multiplied by their emission value)
for pi_intial_prob_pos in 1:num_haplotypes^2
    T1[pi_intial_prob_pos][1] = 1 *  get_emission_value(OBS_ALLELE,pi_intial_prob_pos)
    T2[pi_intial_prob_pos][1] = 0
end

for observation_number in (obs_index + 1):num_SNPs

    HAPS_LINE = split(haplotypes[observation_number])
    OBS_ALLELE = (split(observations[observation_number]))[2]

    for haplotype_state in 1:num_haplotypes^2

        maximum_probability_path_to_point = 0 #WILL THIS WORK?:max(T1[:][observation_number - 2] * ... )
        latest_hap_of_best_path_so_far = 0 #argmax( of above )

        for previous_state in 1:num_haplotypes

            probability_for_state = T1[previous_state][observation_number - 2]*get_transition_value(previous_state,haplotype_state)*get_emission_value(OBS_ALLELE,haplotype_state) #observation_number - 2 because you have -1 since obs starts at 2 and also -1 since we're checking the previous value

            if probability_for_state > maximum_probability_path_to_point

                maximum_probability_path_to_point = probability_for_state
                latest_hap_of_best_path_so_far = previous_state #IS THIS RIGHT - ToDo Check
            end
        end

        T1[haplotype_state][observation_number - 1] = maximum_probability_path_to_point
        T2[haplotype_state][observation_number - 1] = latest_hap_of_best_path_so_far

    end

end

most_likely_path = Array{Int64}(undef, num_SNPs - 1)

most_likely_path[num_SNPs - 1] = argmax(T1[:][num_SNPs - 1])

num_states_remaining = num_SNPs - 1

while num_states_remaining > 1

    most_likely_path[num_states_remaining - 1] = T2[most_likely_path[num_states_remaining]]

end

println(most_likely_path)



#using the convention that the T1 and T2 tables are being filled in with the first label column being AA corresponding to 1,AB corresponding to 2,AC corresponding to 3,AD,...,BA,BB,BC,... etc. meaning that the index can be translated to a unique combo of haplotypes. e.g. if there were 4 haplotypes, ABCD, 6 would be BB since 6%4 is 2 and the floor + haps_index bumps it up one making that also 2.
function get_emission_value(observation, state_index)
    state_tuple = (HAPS_LINE[floor(state_index/num_haplotypes) + hap_index], HAPS_LINE[num_haplotypes%(state_index - 1) + hap_index])

    if (observation == state_tuple[1]) && (observation == state_tuple[2])
        return 0.9
    elseif (observation == state_tuple[1]) || (observation == state_tuple[2])
        return 0.45
    else
        return 0.1
    end
end

function get_transition_value(previous_state_index, state_index)
    previous_state_tuple = (floor(previous_state_index/num_haplotypes) + hap_index,num_haplotypes%(previous_state_index - 1) + hap_index)

    state_tuple = (floor(state_index/num_haplotypes) + hap_index,num_haplotypes%(state_index - 1) + hap_index)

    if (previous_state_tuple[1] == state_tuple[1]) && (previous_state_tuple[2] == state_tuple[2])
        return 0.9
    elseif (previous_state_tuple[1] == state_tuple[1]) || (previous_state_tuple[2] == state_tuple[2])
        return 0.1
    else
        return 0
    end
end
