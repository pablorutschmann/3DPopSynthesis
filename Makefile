OUT="build_output"

build:
	cmake -S . -B ${OUT}
	cmake --build ${OUT}
