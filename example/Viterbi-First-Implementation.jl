#!/usr/bin/env julia
#Joseph Outten - viterbi algorithm implementation
#medium 1/10 of number of observations should be about 10min.
#Shuffle the file to get some subset of the
#Need to parse info file talking about how position can changetransition values - file incorporates function of recombination rate at different nucleotide distances (centimorgans) to see how likely recombination between.
#tell it how recombinant the inidv are to get scalar to mult by the recombination rate - use the viterbi

#Optimizations:
#Declare variable types during declaration
#Get some parallelization
#PASS THESE IN INSTEAD OF A BUNCH OF GLOBAL VARIABLES

#MAY NEED TO MAKE IT A BUNCH OF FUNCTIONS INSTEAD TO REALLY SPEED UP RUNTIME AND AVOID GLOBAL VARIABLES

using DelimitedFiles
using LinearAlgebra

function get_emission_value(observation, typed_prediction, IBD, f_error, t_error) #ToDo CLean this up with maybe passing in an array that can be called for the different cases to make it more readable and concise and possibly faster...

    # println("Observation: ", observation)

    if observation == "0"

        if typed_prediction == "00" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == "0." #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == ".1" #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1." #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".." #DONE

            numerator = dot([1-t_error, 0.5, 0.5, t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        else
            print("Wrong haplotype file format.")

        end

    elseif observation == "1" #(i.e. it's the second or 2 (coded as 1))

        if typed_prediction == "00" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == "0." #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == ".1" #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1." #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".." #DONE

            numerator = dot([t_error, 0.5, 0.5, 1-t_error].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        else
            print("Wrong haplotype file format.")

        end

    elseif observation == "."

        if typed_prediction == "00" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[IBD*(1-f_error)+(1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*f_error + (1-IBD)*(f_error^2)])
            return numerator/denominator

        elseif typed_prediction == "01" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(1-f_error)^2, (1-IBD)*f_error^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "10" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error^2), (1-IBD)*(1-f_error)^2, (1-IBD)*(f_error)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "11" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(IBD*f_error + (1-IBD)*(f_error^2)), (1-IBD)*(f_error)*(1-f_error), (1-IBD)*(f_error)*(1-f_error), IBD*(1-f_error) + (1-IBD)*(1-f_error)^2])
            return numerator/denominator

        elseif typed_prediction == ".0" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == "0." #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*(1-founder_allelic_error), (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*f_error]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        elseif typed_prediction == ".1" #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == "1." #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[(1-IBD)*f_error, (1-IBD)*(1-f_error), (1-IBD)*f_error, (1-IBD)*(1-f_error)])
            return numerator/denominator

        elseif typed_prediction == ".." #DONE

            numerator = dot([1,1,1,1].*[.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            denominator = dot([.25+IBD*.25, .25*IBD, .25*IBD, .25+IBD*.25],[1, 1-IBD, 1-IBD, 1]) #USE SIMD TO VECTORIZE THIS
            return numerator/denominator

        else
            print("Wrong haplotype file format.")
            return 0

        end

    else
        print("Wrong observation file format")
        return 0

    end

    # return 0.9

end

function get_transition_value(previous_state_index, state_index,header_names)

    #take the final column and multiply by 0.001 to get the probability of the
    #recombination in that range, then take the cumulative sum of the probabilities in
    #the column and make a linear model to be able to give two positions and get the
    #recombination probability between them by subtracting their y-values given that
    #model (would just be the cumulative sum if thinking about recombination from the
    #start of the chromosome), then use all this to generate the transition probability
    #Probably good idea to plot it to see if its right by plotting and comparing to
    #the graph that Cory provided
    #Do not make it more likely to recombine if haven't seen one in a while because
    #those are independent events and just looking at one window at a first_time

    #This may need to be updated to reflect that there can be two crossover events in one location


    # if state_index%num_haplotypes == 0
    #     cur_second_value = header_names[num_haplotypes + index_offset_of_first_haplotype]
    # else
    #     cur_second_value = header_names[state_index%(num_haplotypes) + index_offset_of_first_haplotype]
    # end

    # global header_names

    cur_second_value = header_names[(state_index - 1)%(num_haplotypes) + index_offset_of_first_haplotype + 1] #need to do "+1 -1" to ensure that the last item (mod of it == 0) is mapped to the last haplotype combo...Julia array starting at 1 kinda necessitates this

    # if previous_state_index%num_haplotypes == 0
    #     prev_second_value = header_names[num_haplotypes + index_offset_of_first_haplotype]
    # elue
    #     prev_second_value = header_names[previous_state_index%(num_haplotypes) + index_offset_of_first_haplotype]
    # endu
    prev_second_value = header_names[(previous_state_index - 1)%(num_haplotypes) + index_offset_of_first_haplotype + 1]

    state_tuple = (header_names[cld(state_index,(num_haplotypes)) + index_offset_of_first_haplotype], cur_second_value)

    previous_state_tuple = (header_names[cld(previous_state_index,(num_haplotypes)) + index_offset_of_first_haplotype], prev_second_value)

    if previous_state_index == state_index
        return 0.999
    elseif (previous_state_tuple[1] == state_tuple[1]) || (previous_state_tuple[2] == state_tuple[2])
        return 0.001
    else
        return 0
    end

end

first_time = time()
observation_file = open(ARGS[1]) #opening the files
haplotypes_file = open(ARGS[2]) #opening the files

observations = readdlm(observation_file, String) #reading in the lines - May be too big to deal with
test_observations = readdlm(open("alt_observations.dat"), Char)
# println("Beginning of the test file: ", test_observations[1:10])
num_rows = length(observations[:,1])
# positions = observations[:,1]
observations = observations[:,2]
# observations = test_observations
haplotypes = readdlm(haplotypes_file, String) #reading in the readlines
OBS_ALLELE = observations[2] #the observation being dealt with
close(observation_file)
num_SNPs = min(length(observations), length(haplotypes[:,1])) - 1

#PASS THESE IN INSTEAD OF A BUNCH OF GLOBAL VARIABLES
index_offset_of_first_haplotype = 4 #index of first haplotype info minus one (offset)
HAPS_LINE = haplotypes[2,:] #holds the current row to deal with
num_haplotypes = length(HAPS_LINE) - index_offset_of_first_haplotype #holds the number of haplotypes
founder_allelic_error = 0.005
offspring_allelic_error = 0.005


start_of_observations_row = 2 #first row with data
header_names = haplotypes[1,:]

num_hap_combinations = num_haplotypes^2

emission_value_dict = Dict{String,Float64}() #will store the translation from observed (untyped) to probability of being the true genotype. It is of the form emission_value_dict["<observed untyped genotype"] = emission value

most_likely_paths_so_far = Array{Float64}(undef, num_hap_combinations, num_SNPs) #holds the probability of the most likely path so far to this haplotype
most_likely_previous_paths = Array{Float64}(undef, num_hap_combinations, num_SNPs) #holds where the most likely path thus far came from
#
for initial_prob_position in 1:(num_haplotypes^2)

    founder_1 = HAPS_LINE[cld(initial_prob_position, (num_haplotypes)) + index_offset_of_first_haplotype]
    founder_2 = HAPS_LINE[(initial_prob_position - 1)%num_haplotypes + index_offset_of_first_haplotype + 1]
    if founder_1 == founder_2
        IBD = 1
    else
        IBD = 0
    end
    # println("log of the emission: ",  log(100, get_emission_value(string(OBS_ALLELE), string(founder_1, founder_2), IBD, founder_allelic_error, offspring_allelic_error)))
    most_likely_paths_so_far[initial_prob_position,1] = 1 * log(get_emission_value(OBS_ALLELE, founder_1*founder_2, IBD, founder_allelic_error, offspring_allelic_error))
    most_likely_previous_paths[initial_prob_position,1] = 0

end

# print("[")

#findmax(a[:,5])
to_mod = ceil(Int, num_rows/1000)
counter = 0.0

for observation_number in (start_of_observations_row + 1):(num_SNPs + 1)
    # print("num___rows: ", num_rows)
    if (observation_number % to_mod == 0)
        global counter += 1.0
        print("The code is ", counter/10, "% done. ", '\r')
    end
    # if observation_number % 50 == 0
    #     print(".")
    # end

    OBS_ALLELE = observations[observation_number]
    # global
    prev_HAPS_LINE = haplotypes[observation_number - 1,:]
    # global
    HAPS_LINE = haplotypes[observation_number,:]

    for haplotype_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED
        finding_the_way = Array{Float64}(undef, num_hap_combinations)
        for previous_state in 1:(num_haplotypes^2) #OPTIMIZATION NEEDED

            finding_the_way[previous_state] = most_likely_paths_so_far[previous_state,observation_number - 2]+log(get_transition_value(previous_state,haplotype_state, header_names)) #observation_number - 2 because you have -1 since obs starts at 2 and also -1 since we're checking the previous value
            # @time get_transition_value(previous_state, haplotype_state)

        end
        maxval, maxpos = findmax(finding_the_way)
        most_likely_previous_paths[haplotype_state, observation_number - 1] = floor(Int, maxpos)


        founder_1 = HAPS_LINE[cld(haplotype_state, (num_haplotypes)) + index_offset_of_first_haplotype]
        founder_2 = HAPS_LINE[(haplotype_state - 1)%num_haplotypes + index_offset_of_first_haplotype + 1]
        if founder_1 == founder_2
            IBD = 1
        else
            IBD = 0
        end
        most_likely_paths_so_far[haplotype_state, observation_number - 1] = maxval+log(get_emission_value(OBS_ALLELE, founder_1*founder_2, IBD, founder_allelic_error, offspring_allelic_error))
        # @time get_emission_value(Int,OBS_ALLELE, string(founder_1*founder_2), IBD, founder_allelic_error, offspring_allelic_error)

    end

end

# println("]")

prob_max_path, final_hap = findmax(most_likely_paths_so_far[:,num_SNPs])
final_path = Array{String}(undef, num_SNPs)

back_trace = final_hap
back_trace_second = back_trace%num_haplotypes
if back_trace_second == 0
    back_trace_second = num_haplotypes
end
final_path[num_SNPs] = string(header_names[cld(back_trace,num_haplotypes) + index_offset_of_first_haplotype] , " : ", header_names[back_trace_second + index_offset_of_first_haplotype])
for back_trace_step in 1:(num_SNPs - 1)
    global back_trace = floor(Int, most_likely_previous_paths[back_trace, num_SNPs - back_trace_step + 1])

    back_trace_second = back_trace%num_haplotypes

    if back_trace_second == 0
        back_trace_second = num_haplotypes
    end
    final_path[num_SNPs - back_trace_step] = string(header_names[cld(back_trace,num_haplotypes) + index_offset_of_first_haplotype], " : " , header_names[back_trace_second + index_offset_of_first_haplotype])

end

second_time = time()
println("")
println("time elapsed = ", (second_time - first_time))
#println("The end: ", final_path)
haplotypes = haplotypes[:,2]
# println(haplotypes[2], " : " , final_path[1])
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
    # if indexing >= length(final_path)
    #     break
    # end
    # println(haplotypes[indexing+1], " : ", final_path[indexing])
    push!(output, split(final_path[indexing], " : "))
    push!(output[output_indexing], string(haplotypes[indexing+1]))
    global output_indexing += 1
    global prev_final_path = final_path[indexing]
end
# second_time = time()
# println("time elapsed = ", (second_time - first_time))
println("max prob: ", prob_max_path)

println("haplotype", '\t', "start", '\t', "stop", '\t', "ind", '\t', "sex", '\t', "lineID")
prev_entry = ""
line = ""
for haplotype in 1:2

    # global
    line = string(haplotype, '\t', '\t',"1", '\t')
    # global
    prev_entry = output[1][haplotype]
    for row_num in 2:length(output)
        # global output
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
