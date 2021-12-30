all: clean install-report install-protqc
	@echo "Compile the quartet-protqc-report...."
	@bin/lein uberjar
	@printf "\n\n\e[1;32mRun the command for more details: \nsource .env/bin/activate\njava -jar target/uberjar/quartet-protqc-report-*-standalone.jar -h\e[0m"

clean:
	@echo "Clean the environment..."
	@bin/lein clean
	@rm -rf .env .lsp .clj-kondo report/dist report/quartet-proteome-report.egg-info protqc.tar.gz resources/renv/library resources/renv/staging

make-env:
	virtualenv -p python3 .env
	cp resources/bin/protqc.sh .env/bin

install-report: make-env
	. .env/bin/activate
	cd report && python3 setup.py sdist && pip3 install dist/*.tar.gz

install-protqc:
	tar czvf protqc.tar.gz protqc
	R CMD INSTALL protqc.tar.gz
