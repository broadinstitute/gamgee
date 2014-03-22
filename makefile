lib/b2/b2:
	@cd lib/b2 && ./bootstrap.sh && ./b2 install --prefix=../../boost-build && cd ../../

lib/htslib/libhts.a:
	@cd lib/htslib && make lib-static && cd ../../

travis: lib/b2/b2 lib/htslib/libhts.a
	@echo "built all dependencies"

htslib: lib/htslib/libhts.a
	@echo "built htslib"

clean:
	@rm -rf boost-build && cd lib/htslib && make clean && rm ../b2/b2 ../b2/bjam && cd ../../

.PHONY: travis
