1: 01-calc_vp_halo_llist_pbc.o catalogue.o functions.o
	g++ 01-calc_vp_halo_llist_pbc.o catalogue.o functions.o -o 01-calc_vp_halo_llist_pbc.exe -lm -fopenmp
01-calc_vp_halo_llist_pbc.o: 01-calc_vp_halo_llist_pbc.cpp
	g++ -c 01-calc_vp_halo_llist_pbc.cpp -lm -fopenmp

2: 02-calc_vp_halo_llist_pbc_final-var.o catalogue.o functions.o
	g++ 02-calc_vp_halo_llist_pbc_final-var.o catalogue.o functions.o -o 02-calc_vp_halo_llist_pbc_final-var.exe -lm -fopenmp
02-calc_vp_halo_llist_pbc_final-var.o: 02-calc_vp_halo_llist_pbc_final-var.cpp
	g++ -c 02-calc_vp_halo_llist_pbc_final-var.cpp -lm -fopenmp

3: 03-calc_vp_halo_llist_pbc_dyn-var.o catalogue.o functions.o
	g++ 03-calc_vp_halo_llist_pbc_dyn-var.o catalogue.o functions.o -o 03-calc_vp_halo_llist_pbc_dyn-var.exe -lm -fopenmp
03-calc_vp_halo_llist_pbc_dyn-var.o: 03-calc_vp_halo_llist_pbc_dyn-var.cpp
	g++ -c 03-calc_vp_halo_llist_pbc_dyn-var.cpp -lm -fopenmp

4: 04-calc_vp_part_llist_pbc_var.o catalogue.o functions.o
	g++ 04-calc_vp_part_llist_pbc_var.o catalogue.o functions.o -o 04-calc_vp_part_llist_pbc_var.exe -lm -fopenmp
04-calc_vp_part_llist_pbc_var.o: 04-calc_vp_part_llist_pbc_var.cpp
	g++ -c 04-calc_vp_part_llist_pbc_var.cpp -lm -fopenmp


catalogue.o: classes/catalogue.cpp
	g++ -c classes/catalogue.cpp -Wall -Werror -lm
functions.o: classes/functions.cpp
	g++ -c classes/functions.cpp -Wall -Werror -lm
