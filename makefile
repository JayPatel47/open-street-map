build:
	rm -f application.exe
	g++ -std=c++11 -Wall application.cpp dist.cpp osm.cpp tinyxml2.cpp -o application.exe

run:
	./application.exe

buildtest:
	rm -f testing.exe
	g++ -std=c++11 -Wall testing.cpp -o testing.exe

runtest:
	./testing.exe
	
buildtests:
	rm -f tests.exe
	g++ -std=c++11 -Wall tests.cpp -o tests.exe
	
runtests:
	./tests.exe	

valgrind:
	valgrind --tool=memcheck --leak-check=yes ./application.exe
