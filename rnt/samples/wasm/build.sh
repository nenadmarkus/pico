cat ../../cascades/facefinder | hexdump -v -e '16/1 "0x%x," "\n"' > facefinder.hex
emcc main.c -o wasmpico.js -O3 -s EXPORTED_FUNCTIONS="['_find_faces', '_malloc', '_free']" -s WASM=1
rm facefinder.hex
