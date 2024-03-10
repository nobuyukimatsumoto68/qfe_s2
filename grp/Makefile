all: gen_grp_o3

.phony: gen_grp_o3 clean

clean:
	$(RM) -rf elem* lattice* site_g*
	rm -rf __pycache__

gen_grp_o3: gen_grp_o3.py grp_o3.py quat.py
	python3 gen_grp_o3.py
