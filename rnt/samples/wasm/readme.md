# Compiling `pico` to WebAssembly

First, run

	bash build.sh

to generate the Wasm binary an associated JS boilerplate.

Next, download the webcam library from <https://github.com/nenadmarkus/picojs> and place it in this directory:

	wget https://raw.githubusercontent.com/nenadmarkus/picojs/master/examples/camvas.js

Now you can open `index.html` in your web browser and try how it works in practice.
Note that some browsers won't allow running the demo locally (without a webserver).

Related links:

* post explaining the details: <https://nenadmarkus.com/p/pico-to-wasm/>;
* webcam demo: <https://nenadmarkus.com/p/pico-to-wasm/demo/>.
