PACKET_SIZE="15000"
TRACE_PATH="../data/trace_file/synthetic/multi_process/"

#----------------------------------------------Fermi-------------------------------------------------

make clean; 
make KT=Fermi_kernel;

echo "Fermi_2048patterns_kernel"

STRUCTURAL_CONFIG="--automata_file ../data/mem_layout/Fermi_2048patterns_mem.txt"

OCCUPANCY_STREAM_NUM="75" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '

OUTFILE='../data/out_file/Fermi_traversal_0.50prob.out'
TRACE_CONFIG='--pkt_num 1 --trace_num 8 --tracefile_names '$TRACE_PATH'Fermi_2048patterns_depth_s0_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s1_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s2_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s3_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s4_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s5_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s6_p0.50.trace '$TRACE_PATH'Fermi_2048patterns_depth_s7_p0.50.trace'

echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1


#----------------------------------------------SPM-------------------------------------------------


make clean; 
make KT=SPM_kernel;

echo "SPM_2048patterns_kernel"
STRUCTURAL_CONFIG="--automata_file ../data/mem_layout/SPM_2048patterns_mem.txt"

OCCUPANCY_STREAM_NUM="60" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '



OUTFILE='../data/out_file/SPM_traversal_0.50prob.out'
TRACE_CONFIG='--pkt_num 1 --trace_num 8 --tracefile_names '$TRACE_PATH'SPM_2048patterns_depth_s0_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s1_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s2_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s3_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s4_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s5_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s6_p0.50.trace '$TRACE_PATH'SPM_2048patterns_depth_s7_p0.50.trace'
	echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1"
	./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1



#----------------------------------------------ER-------------------------------------------------


make clean; 
make KT=ER_kernel;

echo "ER_2048patterns_kernel"

STRUCTURAL_CONFIG="--automata_file ../data/mem_layout/ER_2048patterns_mem.txt"

OCCUPANCY_STREAM_NUM="15" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '


OUTFILE='../data/out_file/ER_traversal_0.50prob.out'
TRACE_CONFIG='--pkt_num 1 --trace_num 8 --tracefile_names '$TRACE_PATH'ER_2048patterns_depth_s0_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s1_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s2_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s3_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s4_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s5_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s6_p0.50.trace '$TRACE_PATH'ER_2048patterns_depth_s7_p0.50.trace'

echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE 2>&1


#---------------------------------------------------------------------MAKING NFAS-------------------------------------------------------------------
PATTERN_CONFIG="2048"
TRACE_CONFIG="--trace_num 1 --tracefile_names ../data/trace_file/stream_Ecym_9184_1024ups.tsv --pkt_num 1"

#----------------------------------------------HAmming-------------------------------------------------
make clean; 
make KT=Hamming_k8_d1_kernel IT=DIVIDE_TRACE;

echo 'Hamming_k8_d1_'$PATTERN_CONFIG'patterns'
OUTFILE="../data/out_file/hamming_k8_d1_traversal.out"
STRUCTURAL_CONFIG='--automata_file ../data/mem_layout/Hamming_k8_d1_'$PATTERN_CONFIG'patterns_mem.txt'

OCCUPANCY_STREAM_NUM="82" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '

echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE


make clean; 
make KT=Hamming_k20_d3_kernel IT=DIVIDE_TRACE;

echo 'Hamming_k20_d3_'$PATTERN_CONFIG'patterns'
OUTFILE="../data/out_file/hamming_k20_d3_traversal.out"
STRUCTURAL_CONFIG='--automata_file ../data/mem_layout/Hamming_k20_d3_'$PATTERN_CONFIG'patterns_mem.txt'

OCCUPANCY_STREAM_NUM="15" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '

echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE

#----------------------------------------------Levenshtein-------------------------------------------------
make clean; 
make KT=Levenshtein_k8_d1_kernel IT=DIVIDE_TRACE;

echo 'levenshtein_k8_d1_'$PATTERN_CONFIG'patterns'
OUTFILE="../data/out_file/levenshtein_k8_d1_traversal.out"
STRUCTURAL_CONFIG='--automata_file ../data/mem_layout/Levenshtein_k8_d1_'$PATTERN_CONFIG'patterns_mem.txt'

OCCUPANCY_STREAM_NUM="75" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '


echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG> $OUTFILE"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE # -k 8 -dist 1



make clean; 
make KT=Levenshtein_k20_d3_kernel IT=DIVIDE_TRACE;

echo 'levenshtein_k20_d3_'$PATTERN_CONFIG'patterns'
OUTFILE="../data/out_file/levenshtein_k20_d3_traversal.out"
STRUCTURAL_CONFIG='--automata_file ../data/mem_layout/Levenshtein_k20_d3_'$PATTERN_CONFIG'patterns_mem.txt'

OCCUPANCY_STREAM_NUM="15" #in reality its the number of streams
OCCUPANCY_STREAM_NUM_CONFIG='--stream_count '$OCCUPANCY_STREAM_NUM' '

echo "./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE"
./regex_gpu $STRUCTURAL_CONFIG $TRACE_CONFIG --pkt_size $PACKET_SIZE --device 0 $OCCUPANCY_STREAM_NUM_CONFIG > $OUTFILE
