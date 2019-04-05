#Joseph Outten - viterbi algorithm implementation

function get_emission_value(observation, state_index)

    global HAPS_LINE
    #println("Heres a HAPS_LINE in emission: ", HAPS_LINE)
    if state_index%num_haplotypes == 0
        second_value = HAPS_LINE[num_haplotypes + to_add_to_first_haplotype]
    else
        second_value = HAPS_LINE[state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    state_tuple = (HAPS_LINE[cld(state_index,(num_haplotypes)) + to_add_to_first_haplotype], second_value)
    #println("Emission State Tuple")
    #println(string(state_index, " : ", state_tuple))
    return floor(Int, state_index)

end

function get_transition_value(previous_state_index, state_index)

    global HAPS_LINE
    global prev_HAPS_LINE

    if state_index%num_haplotypes == 0
        cur_second_value = HAPS_LINE[num_haplotypes + to_add_to_first_haplotype]
    else
        cur_second_value = HAPS_LINE[state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    if previous_state_index%num_haplotypes == 0
        prev_second_value = prev_HAPS_LINE[num_haplotypes + to_add_to_first_haplotype]
    else
        prev_second_value = prev_HAPS_LINE[previous_state_index%(num_haplotypes) + to_add_to_first_haplotype]
    end

    state_tuple = (HAPS_LINE[cld(state_index,(num_haplotypes)) + to_add_to_first_haplotype], cur_second_value)

    previous_state_tuple = (prev_HAPS_LINE[cld(previous_state_index,(num_haplotypes)) + to_add_to_first_haplotype], prev_second_value)
    #println("prev_HAPS_LINE: ", prev_HAPS_LINE)
    #println("HAPS_LINE : ", HAPS_LINE)
    #print("Transition State Tuple ")
    #print(string(state_index, " : ", state_tuple, " , "))
    #println(string(previous_state_index, " : ", previous_state_tuple, " end "))

    return floor(Int, previous_state_index)
end

#Need to make the prev_HAPS_LINE and update both of them during the algorithm
first_time = time()
observation_file = open("observations.dat") #opening the files
haplotypes_file = open("test.dat") #opening the files

observations = readlines(observation_file) #reading in the lines - May be too big to deal with
haplotypes = readlines(haplotypes_file) #reading in the readlines
OBS_ALLELE = (split(observations[2]))[2] #the observation being dealt with

global to_add_to_first_haplotype = 4 #index of first haplotype info minus one (offset)
global HAPS_LINE = split(haplotypes[2]) #holds the current row to deal with
global prev_HAPS_LINE = split(haplotypes[2]) #holds the previous row to deal with
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
#(15, 3)

for observation_number in (start_of_observations_row + 1):(num_SNPs + 1)

    #println("num_SNPs: ", num_SNPs)
    #println("current row: ", observation_number)

    #if observation_number % floor(Int, num_SNPs/1) == 0
        print(".")
    #end

    OBS_ALLELE = (split(observations[observation_number]))[2]
    global prev_HAPS_LINE = split(haplotypes[observation_number - 1])
    global HAPS_LINE = split(haplotypes[observation_number])

    #println("Updated prev_HAPS_LINE: ", prev_HAPS_LINE)
    #println("Updated HAPS_LINE: ", HAPS_LINE)

    for haplotype_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED
        finding_the_way = Array{Float64}(undef, (num_haplotypes)^2)
        for previous_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED

            finding_the_way[previous_state] = most_likely_paths_so_far[previous_state,observation_number - 2]*get_transition_value(previous_state,haplotype_state) #observation_number - 2 because you have -1 since obs starts at 2 and also -1 since we're checking the previous value

        end

        #println("finding the way outside: ", finding_the_way)
        maxval, maxpos = findmax(finding_the_way)
        most_likely_previous_paths[haplotype_state, observation_number - 1] = floor(Int, maxpos)
        most_likely_paths_so_far[haplotype_state, observation_number - 1] = maxval*get_emission_value(OBS_ALLELE, haplotype_state)

    end

end

println("]")

#println("most_likely_previous_paths: ", most_likely_previous_paths)
#println("most_likely_paths_so_far: ", most_likely_paths_so_far)

prob_max_path, final_hap = findmax(most_likely_paths_so_far[:,num_SNPs])
final_path = Array{Int64}(undef, num_SNPs)

back_trace = final_hap
final_path[num_SNPs] = back_trace
for back_trace_step in 1:(num_SNPs - 1)
    #println("back_trace: ", back_trace)
    global back_trace = floor(Int, most_likely_previous_paths[back_trace, num_SNPs - back_trace_step + 1])
    final_path[num_SNPs - back_trace_step] = back_trace

end

#println("The end: ", final_path)
second_time = time()
println("time elapsed = ", (second_time - first_time))
