#Joseph Outten - viterbi algorithm implementation
#medium 1/10 of number of observations should be about 10min.
#Shuffle the file to get some subset of the
#Need to parse info file talking about how position can changetransition values - file incorporates function of recombination rate at different nucleotide distances (centimorgans) to see how likely recombination between.
#tell it how recombinant the inidv are to get scalar to mult by the recombination rate - use the viterbi

function get_emission_value(observation, state_index)

    global HAPS_LINE
    if state_index%num_haplotypes == 0
        second_value = HAPS_LINE[num_haplotypes + to_add_to_first_haplotype]
    else
        second_value = HAPS_LINE[state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    temp_state_tuple = (HAPS_LINE[cld(state_index,(num_haplotypes)) + to_add_to_first_haplotype], second_value)

    if temp_state_tuple[1] == "."
        temp_1 = "0"
    else
        temp_1 = temp_state_tuple[1]

    end
    if temp_state_tuple[2] == "."
        temp_2 = "0"
    else
        temp_2 = temp_state_tuple[2]
    end
    state_tuple = (temp_1, temp_2)
    return floor(Int, parse(Int,state_tuple[1]) + parse(Int, state_tuple[2]))

end

function get_transition_value(previous_state_index, state_index)

    global header_names

    if state_index%num_haplotypes == 0
        cur_second_value = header_names[num_haplotypes + to_add_to_first_haplotype]
    else
        cur_second_value = header_names[state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    if previous_state_index%num_haplotypes == 0
        prev_second_value = header_names[num_haplotypes + to_add_to_first_haplotype]
    else
        prev_second_value = header_names[previous_state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    state_tuple = (header_names[cld(state_index,(num_haplotypes)) + to_add_to_first_haplotype], cur_second_value)

    previous_state_tuple = (header_names[cld(previous_state_index,(num_haplotypes)) + to_add_to_first_haplotype], prev_second_value)

    if previous_state_index == state_index
        return 0.9
    else
        return 0.01
    end
end

first_time = time()
observation_file = open(ARGS[1]) #opening the files
haplotypes_file = open(ARGS[2]) #opening the files

observations = readlines(observation_file) #reading in the lines - May be too big to deal with
haplotypes = readlines(haplotypes_file) #reading in the readlines
OBS_ALLELE = (split(observations[2]))[2] #the observation being dealt with

global to_add_to_first_haplotype = 4 #index of first haplotype info minus one (offset)
global HAPS_LINE = split(haplotypes[2]) #holds the current row to deal with
global num_haplotypes = length(HAPS_LINE) - to_add_to_first_haplotype #holds the number of haplotypes

start_of_observations_row = 2 #first row with data
num_SNPs = length(haplotypes) - 1
header_names = split(haplotypes[1])

most_likely_paths_so_far = Array{Float64}(undef, (num_haplotypes)^2, num_SNPs) #holds the probability of the most likely path so far to this haplotype
most_likely_previous_paths = Array{Float64}(undef, (num_haplotypes)^2, num_SNPs) #holds where the most likely path thus far came from

for initial_prob_position in 1:(num_haplotypes^2)

    most_likely_paths_so_far[initial_prob_position,1] = 1 * get_emission_value(OBS_ALLELE, initial_prob_position)
    most_likely_previous_paths[initial_prob_position,1] = 0

end

print("[")

#findmax(a[:,5])

for observation_number in (start_of_observations_row + 1):(num_SNPs + 1)

    if observation_number % floor(Int, num_SNPs/50) == 0
        print(".")
    end

    OBS_ALLELE = (split(observations[observation_number]))[2]
    global prev_HAPS_LINE = split(haplotypes[observation_number - 1])
    global HAPS_LINE = split(haplotypes[observation_number])

    for haplotype_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED
        finding_the_way = Array{Float64}(undef, (num_haplotypes)^2)
        for previous_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED

            finding_the_way[previous_state] = most_likely_paths_so_far[previous_state,observation_number - 2]*get_transition_value(previous_state,haplotype_state) #observation_number - 2 because you have -1 since obs starts at 2 and also -1 since we're checking the previous value

        end

        maxval, maxpos = findmax(finding_the_way)
        most_likely_previous_paths[haplotype_state, observation_number - 1] = floor(Int, maxpos)
        most_likely_paths_so_far[haplotype_state, observation_number - 1] = maxval*get_emission_value(OBS_ALLELE, haplotype_state)

    end

end

println("]")

prob_max_path, final_hap = findmax(most_likely_paths_so_far[:,num_SNPs])
final_path = Array{String}(undef, num_SNPs)

back_trace = final_hap
back_trace_second = back_trace%num_haplotypes
if back_trace_second == 0
    back_trace_second = num_haplotypes
end
final_path[num_SNPs] = string(header_names[cld(back_trace,num_haplotypes) + to_add_to_first_haplotype] , " : ", header_names[back_trace_second + to_add_to_first_haplotype])
for back_trace_step in 1:(num_SNPs - 1)
    global back_trace = floor(Int, most_likely_previous_paths[back_trace, num_SNPs - back_trace_step + 1])

    back_trace_second = back_trace%num_haplotypes

    if back_trace_second == 0
        back_trace_second = num_haplotypes
    end
    final_path[num_SNPs - back_trace_step] = string(header_names[cld(back_trace,num_haplotypes) + to_add_to_first_haplotype], " : " , header_names[back_trace_second + to_add_to_first_haplotype])

end

#println("The end: ", final_path)
println(split(haplotypes[2])[2], " : " , final_path[1])
indexing = 2
prev_final_path = final_path[1]
while indexing < length(final_path)
    while (final_path[indexing] == prev_final_path)
        global indexing += 1
        if indexing >= length(final_path)
            break
        end
    end
    println(split(haplotypes[indexing+1])[2], " : ", final_path[indexing])
    global prev_final_path = final_path[indexing]
end
second_time = time()
println("time elapsed = ", (second_time - first_time))
println("max prob: ", prob_max_path)
