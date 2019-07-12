#!/usr/bin/env julia
#Joseph Outten - viterbi algorithm implementation

using DelimitedFiles
using LinearAlgebra

function get_emission_value(observation, typed_prediction, IBD, f_error, t_error)

    """
    Returns the emission value for a haplotype state given the observation,
    haplotype state typed prediction, identical by descent, founder typing
    error, and individual typing error. The equations used to generate the
    probabilities were drawn from Tables S1 and S2 from the following paper:

    Zheng C, Boer MP, van Eeuwijk FA. Reconstruction of Genome Ancestry Blocks in Multiparental Populations. Genetics. 2015;200(4):1073–1087. doi:10.1534/genetics.115.177873
    """

    if observation == "0" #i.e. observation is the reference allele

        if typed_prediction == "00" #i.e. typed prediction is homozygous for the reference

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01"

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10"

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11"

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0"

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == "0."

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == ".1"

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1."

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".."

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            return numerator/denominator

        else
            print("Wrong haplotype file format.")

        end

    elseif observation == "1"

        if typed_prediction == "00"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == "0."

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == ".1"

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1."

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".."

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            return numerator/denominator

        else
            print("Wrong haplotype file format.")

        end

    elseif observation == "."

        if typed_prediction == "00"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == "0."

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error])
            return numerator/denominator

        elseif typed_prediction == ".1"

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1."

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".."

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1])
            return numerator/denominator

        else
            print("Wrong haplotype file format.")
            return 0

        end

    else
        print("Wrong observation file format")
        return 0

    end

end

function get_transition_value(previous_state_index, state_index, num_haplotypes, index_offset_of_first_haplotype)

    """
    Returns the probability of one haplotype state transitioning to another
    between sequence reads given the previous state index and current state
    index (labels for those respective states) as well as the number of
    relevant haplotypes and the offset of the first haplotype in the file.
    """

    cur_second_value = (state_index - 1)%(num_haplotypes) + index_offset_of_first_haplotype + 1

    prev_second_value = (previous_state_index - 1)%(num_haplotypes) + index_offset_of_first_haplotype + 1

    state_tuple = (cld(state_index,(num_haplotypes)) + index_offset_of_first_haplotype, cur_second_value)

    previous_state_tuple = (cld(previous_state_index,(num_haplotypes)) + index_offset_of_first_haplotype, prev_second_value)

    if previous_state_index == state_index
        return 0.999
    elseif (previous_state_tuple[1] == state_tuple[1]) || (previous_state_tuple[2] == state_tuple[2])
        return 0.001
    else
        return 0
    end

end

println("Using ", Threads.nthreads(), " threads.")

#log the beginning of the program for timing purposes
first_time = time()

#loading the files given as command line arguments
observation_file = open(ARGS[1]) #First argument: single individual sequencing data to be imputed. Must be tab separated with two columns under a header. The first is a column of positions and the second is a column of observations (0, 1, or . since we assume no heterozygotes).
haplotypes_file = open(ARGS[2]) #Second argument: vcf file of all the founder haplotypes.

#setting up necessary parameters
observations = readdlm(observation_file, String)[:,2]
haplotypes = readdlm(haplotypes_file, String)
close(observation_file)
close(haplotypes_file)
OBS_ALLELE = observations[2]
num_SNPs = min(length(observations), length(haplotypes[:,1])) - 1
index_offset_of_first_haplotype = 4
HAPS_LINE = haplotypes[2,:]
num_haplotypes = length(HAPS_LINE) - index_offset_of_first_haplotype
founder_allelic_error = 0.005
offspring_allelic_error = 0.005
start_of_observations_row = 2
header_names = haplotypes[1,:]
num_hap_combinations = num_haplotypes^2
most_likely_paths_so_far = Array{Float64}(undef, num_hap_combinations, num_SNPs) #holds the probability of the most likely path so far to this haplotype
most_likely_previous_paths = Array{Float64}(undef, num_hap_combinations, num_SNPs) #holds where the most likely path thus far came from

#initializes the first column of most_likely_paths_so_far with the appropriate emission values.
#Threads.@threads
for initial_prob_position in 1:(num_haplotypes^2)

    founder_1 = HAPS_LINE[cld(initial_prob_position, (num_haplotypes)) + index_offset_of_first_haplotype]
    founder_2 = HAPS_LINE[(initial_prob_position - 1)%num_haplotypes + index_offset_of_first_haplotype + 1]
    if founder_1 == founder_2
        IBD = 1
    else
        IBD = 0
    end
    most_likely_paths_so_far[initial_prob_position,1] = log(get_emission_value(OBS_ALLELE, founder_1*founder_2, IBD, founder_allelic_error, offspring_allelic_error))
    most_likely_previous_paths[initial_prob_position,1] = 0

end

#stuff for reporting code completion
to_mod = ceil(Int, num_SNPs/1000)
counter = 0.0

#implements the Viterbi algorithm generating emission and transition values on the fly
for observation_number in (start_of_observations_row + 1):(num_SNPs + 1)

    #prints the approximate completion status of the code
    if (observation_number % to_mod == 0)

        global counter += 1.0
        print("The code is ", counter/10, "% done. ", '\r')

    end

    OBS_ALLELE = observations[observation_number]
    HAPS_LINE = haplotypes[observation_number,:]

#Threads.@threads
    for haplotype_state in 1:(num_haplotypes^2)

        finding_the_way = Array{Float64}(undef, num_hap_combinations)

        for previous_state in 1:(num_haplotypes^2)

            finding_the_way[previous_state] = most_likely_paths_so_far[previous_state,observation_number - 2]+log(get_transition_value(previous_state,haplotype_state, num_haplotypes, index_offset_of_first_haplotype))

        end

        maxval, maxpos = findmax(finding_the_way)
        most_likely_previous_paths[haplotype_state, observation_number - 1] = floor(Int, maxpos)
        founder_1 = cld(haplotype_state, (num_haplotypes)) + index_offset_of_first_haplotype
        founder_2 = (haplotype_state - 1)%num_haplotypes + index_offset_of_first_haplotype + 1

        if founder_1 == founder_2
            IBD = 1
        else
            IBD = 0
        end

        most_likely_paths_so_far[haplotype_state, observation_number - 1] = maxval+log(get_emission_value(OBS_ALLELE, HAPS_LINE[founder_1]*HAPS_LINE[founder_2], IBD, founder_allelic_error, offspring_allelic_error))

    end

end

#backtracking through the most_likely_previous_paths to find the highest likelihood path
prob_max_path, final_hap = findmax(most_likely_paths_so_far[:,num_SNPs])
final_path = Array{String}(undef, num_SNPs)
back_trace = final_hap
back_trace_second = back_trace%num_haplotypes

if back_trace_second == 0
    back_trace_second = num_haplotypes
end

#formatting final_path to hold the name of the appropriate haplotype at each position
final_path[num_SNPs] = string(header_names[cld(back_trace,num_haplotypes) + index_offset_of_first_haplotype] , " : ", header_names[back_trace_second + index_offset_of_first_haplotype])

for back_trace_step in 1:(num_SNPs - 1)

    global back_trace = floor(Int, most_likely_previous_paths[back_trace, num_SNPs - back_trace_step + 1])
    back_trace_second = back_trace%num_haplotypes

    if back_trace_second == 0
        back_trace_second = num_haplotypes
    end

    final_path[num_SNPs - back_trace_step] = string(header_names[cld(back_trace,num_haplotypes) + index_offset_of_first_haplotype], " : " , header_names[back_trace_second + index_offset_of_first_haplotype])

end

#printing time elapsed
second_time = time()
println("\ntime elapsed = ", (second_time - first_time))

#formatting final output
haplotypes = haplotypes[:,2]
output = Array{String}[]
push!(output, split(final_path[1], " : "))
push!(output[1], "1")
output_indexing = 2
indexing = 2
prev_final_path = final_path[1]

while indexing < length(final_path)

    while (final_path[indexing] == prev_final_path)

        global indexing += 1
        if indexing >= length(final_path)
            break
        end

    end

    push!(output, split(final_path[indexing], " : "))
    push!(output[output_indexing], string(haplotypes[indexing+1]))
    global output_indexing += 1
    global prev_final_path = final_path[indexing]

end

#printing final output
println("max prob: ", prob_max_path)

println("haplotype", '\t', "start", '\t', "stop", '\t', "ind", '\t', "sex", '\t', "lineID")
prev_entry = ""
line = ""
for haplotype in 1:2

    line = string(haplotype, '\t', '\t',"1", '\t')
    prev_entry = output[1][haplotype]

    for row_num in 2:length(output)

        entry = output[row_num][haplotype]
        if entry != prev_entry
            line = line * string(output[row_num][3], "   ", '\t', "-", '\t', "-", '\t', prev_entry)
            println(line)
            line = string(haplotype, "  ",  '\t', "   " ,output[row_num][3], '\t')
            prev_entry = entry
        end
        
    end
    line = line * string(output[length(output)][3], "   ", '\t', "-", '\t', "-", '\t', prev_entry)
    println(line)
end
