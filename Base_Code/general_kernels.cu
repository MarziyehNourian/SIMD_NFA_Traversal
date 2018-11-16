//By: Marziyeh Nourian
#include "general_kernels.h"

__global__ void fixed_topology_kernel(  	unsigned* result_bit_vector,
											char* stream_sequences,
											#ifdef RODC_ON
											struct preprocessed_full_reference_char_sequence const * __restrict__  preprocessed_input,
											#else
											struct preprocessed_full_reference_char_sequence * preprocessed_input,
											#endif
											unsigned bit_chunks_per_state_vector, //number of integers in sv
											unsigned char_filled_ints_per_packet,
											unsigned num_packets,
											unsigned warp_efficient_stream_count,
											unsigned occupancy_efficient_stream_count,
											unsigned ref_block_count,
											unsigned batch_count,	//the number of word-chunks (32-bit) that holds the symbols of NFA s on each block
											unsigned accepting_states_count  ){

	//Registers
	unsigned ctof =							bit_chunks_per_state_vector;
	unsigned * u_stream_sequences = 	  			(unsigned *)stream_sequences;
	unsigned packet_idx = 						0;
	unsigned p_idx =						0;	
	unsigned char_filled_int_idx =					0;	
	char 	s_char =						0;	
	unsigned s_chars = 						0;	
	int s_char_index = 						0;

	//id of the input stream based on the thread and block count
	int stream_id = warp_efficient_stream_count*(blockIdx.x/ref_block_count) + (threadIdx.x/batch_count);
	if (stream_id > (occupancy_efficient_stream_count*warp_efficient_stream_count)-1 ) printf("stream-id exceed the allowed range!\n");

	extern __shared__ unsigned cfBV[];

	for(packet_idx=0; packet_idx < char_filled_ints_per_packet*num_packets; packet_idx+=char_filled_ints_per_packet){
	//set first state (i.e. non-anchored start state)
    	cfBV[threadIdx.x] = 0xFFFFFFFF;
	    cfBV[ctof + threadIdx.x] = 0xFFFFFFFF;

		for(char_filled_int_idx = packet_idx; char_filled_int_idx < packet_idx + char_filled_ints_per_packet; char_filled_int_idx += 1){

			s_chars = u_stream_sequences[(stream_id*char_filled_ints_per_packet*num_packets)+char_filled_int_idx];

			for(s_char_index = 3; s_char_index >= 0; s_char_index -= 1){

				s_char = (char)(  (s_chars>>(8*(3-s_char_index))) & 0x000000FFu   );

				topology_specific_traversal(s_char,preprocessed_input,cfBV,ctof,bit_chunks_per_state_vector,accepting_states_count,ref_block_count,batch_count);

			}

		}

		//Write to global memory
		fill_results_array(result_bit_vector, cfBV, ctof, bit_chunks_per_state_vector, accepting_states_count , packet_idx, &p_idx);

	}
}


__device__ unsigned match_check(char s_char,
								unsigned passed_chars, //is calculated by the compiler and depends on the src and dst of the TX
								unsigned tx_char_count,
								#ifdef RODC_ON
								struct preprocessed_full_reference_char_sequence const * __restrict__ preprocessed_input,
								#else
								struct preprocessed_full_reference_char_sequence * preprocessed_input,
								#endif
								unsigned ref_block_count,
								unsigned batch_count){//n is the number of 32-bit words that contain all N nfas
		unsigned mask_res  = 	0x00000000;

		int tx_char_iterator;
		for(tx_char_iterator = 0; tx_char_iterator < tx_char_count; tx_char_iterator++){

		unsigned mask  = 	0x00000000;
		unsigned offset_in_SOA = passed_chars*batch_count;

		unsigned p_chars_7 = 	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_7_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_6 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_6_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_5 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_5_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_4 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_4_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_3 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_3_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_2 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_2_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_1 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_1_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];
		unsigned p_chars_0 =	preprocessed_input[(blockIdx.x%ref_block_count)].SOA_chunks_0_word[offset_in_SOA+(batch_count*tx_char_iterator)+(threadIdx.x%batch_count)];

		int j = 0;
		char p_char;

		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_7>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_6>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_5>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_4>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_3>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_2>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_1>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask = mask<<4;
		for(j=3; j>=0; j-=1){
			p_char = (char)((p_chars_0>>(8*(3-j))) & 0x000000FFu);
			if(s_char==p_char){
				mask |= ( (0x1) << j );
			}
		}
		mask_res |= mask;
	}
	return mask_res;
}

__device__ void character_transitions_update(unsigned * cfBV, unsigned ctof,unsigned src, unsigned dst, unsigned mask,unsigned has_wildcard,unsigned has_neg,unsigned has_positive){

	unsigned current_index = (blockDim.x*src) + threadIdx.x;
	unsigned future_index = (blockDim.x*dst) + threadIdx.x;

	unsigned current = cfBV[current_index];

	if(has_wildcard){
		cfBV[ctof + future_index] |= current;
	}else if(has_positive){
		cfBV[ctof + future_index] |= mask & current;
	}else if(has_neg){
		cfBV[ctof + future_index] |= ~(mask) & current;
	}
}

__device__ void update_StateVector(unsigned * cfBV, unsigned ctof, unsigned bit_chunks_per_state_vector, unsigned accepting_states_count){
	int i=0;
	for(i=threadIdx.x + blockDim.x; i<bit_chunks_per_state_vector-((accepting_states_count)*blockDim.x); i+=blockDim.x){
		cfBV[i] = cfBV[ctof + i];
		cfBV[ctof + i] = 0;
	}//update current with future loop
	for(i=threadIdx.x + (bit_chunks_per_state_vector-((accepting_states_count)*blockDim.x)); i<bit_chunks_per_state_vector; i+=blockDim.x){
		cfBV[i] |= cfBV[ctof + i];
		cfBV[ctof + i] = 0;
	}//accumulative update of accepting states
}

__device__ void fill_results_array(unsigned * result_bit_vector, unsigned * cfBV, unsigned ctof, unsigned bit_chunks_per_state_vector, unsigned accepting_states_count, unsigned packet_idx, unsigned * p_idx){
	int i = 0;
	#ifdef STATE_VECTOR_DEBUG
	if( (packet_idx==0)){
		for(i=threadIdx.x; i<bit_chunks_per_state_vector; i+=blockDim.x){
			result_bit_vector[bit_chunks_per_state_vector*blockIdx.x + i] = cfBV[i];
			cfBV[i]=0;
		}
	}
	#else //Copy only accepting state vector portions


	unsigned offset = (*p_idx)*(accepting_states_count*gridDim.x*blockDim.x) + ( accepting_states_count*blockDim.x*blockIdx.x );
	for(i=(bit_chunks_per_state_vector - (accepting_states_count*blockDim.x) ) + threadIdx.x; i<bit_chunks_per_state_vector; i+=blockDim.x){
		result_bit_vector[ offset + i-((bit_chunks_per_state_vector - (accepting_states_count*blockDim.x) ))] = cfBV[i];
	}

	for(i=threadIdx.x; i<2*bit_chunks_per_state_vector; i+=blockDim.x){
		cfBV[i] = 0; 
	}
	(*p_idx)++;
	#endif
}


