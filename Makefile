.PHONY: clean planned springs view_planned view_springs write_planned write_springs

default: gridify planned springs

# ----------------------------------------

gridify:
	make write_gridify
	make view_gridify

planned:
	make write_planned
	make view_planned

springs:
	make write_springs
	make view_springs

clean:
	rm *.dot *.png

# ----------------------------------------

view_gridify:
	open states_gridify.png

view_planned:
	open states_planned.png

view_springs:
	open states_springs.png

write_gridify:
	./states_viz.py --gridify
	neato -n2 -Tpng -o states_gridify.png states_gridify.dot

write_planned:
	./states_viz.py --planned
	neato -n2 -Tpng -o states_planned.png states_planned.dot

write_springs:
	./states_viz.py --springs
	dot -Kfdp -n1 -Tpng -o states_springs.png states_springs.dot
