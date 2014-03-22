lib/b2/b2:
	@cd lib/b2 && ./bootstrap.sh && cd ../../

lib/htslib/libhts.a:
	@cd lib/htslib && make lib-static && cd ../../

dependencies: lib/b2/b2 lib/htslib/libhts.a
	@echo "built all dependencies"

clean:
	@cd lib/htslib && make clean && rm ../b2/b2 ../b2/bjam && cd ../../

.PHONY: dependencies
