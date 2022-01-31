OUT="build_output"

build:
	mkdir ${OUT}
	cmake -S . -B ${OUT}
	cmake --build ${OUT}
	