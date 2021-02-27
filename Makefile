.PHONY: all clean
.PHONY: nations_gridify nations_planned
.PHONY: states_gridify states_planned states_springs

.PHONY: view_nations_gridify write_nations_gridify
.PHONY: view_nations_planned write_nations_planned

.PHONY: view_states_gridify write_states_gridify
.PHONY: view_states_planned write_states_planned
.PHONY: view_states_springs write_states_springs

default: states_planned

all: states_planned states_springs states_gridify nations_planned nations_gridify

# ----------------------------------------

# Note: Intentionally not removing dot files
clean:
	rm *.dot *.png

nations_gridify:
	make write_nations_gridify
	make view_nations_gridify

nations_planned:
	make write_nations_planned
	make view_nations_planned

states_gridify:
	make write_states_gridify
	make view_states_gridify

states_planned:
	make write_states_planned
	make view_states_planned

states_springs:
	make write_states_springs
	make view_states_springs

# ----------------------------------------

view_nations_gridify:
	open nations_gridify.png

view_nations_planned:
	open nations_planned.png

view_states_gridify:
	open states_gridify.png

view_states_planned:
	open states_planned.png

view_states_springs:
	open states_springs.png

# ----------------------------------------

write_nations_gridify:
	./states_viz.py --gridify --nations --coords_file nations.latlong --dotfile=nations_gridify.dot
	neato -n2 -Tpng -o nations_gridify.png nations_gridify.dot

write_nations_planned:
	./states_viz.py --planned --nations --coords_file nations.latlong --dotfile=nations_planned.dot
	neato -n2 -Tpng -o nations_planned.png nations_planned.dot

write_states_gridify:
	./states_viz.py --gridify  --dotfile=states_gridify.dot
	neato -n2 -Tpng -o states_gridify.png states_gridify.dot

write_states_planned:
	./states_viz.py --planned --coords_file states_v3.coords --dotfile=states_planned.dot
	neato -n2 -Tpng -o states_planned.png states_planned.dot

write_states_springs:
	./states_viz.py --springs --dotfile=states_springs.dot
	dot -Kfdp -Tpng -o states_springs.png states_springs.dot

# ----------------------------------------

test:
	./states_viz.py --test
