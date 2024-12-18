FROM python:latest

RUN \
	git clone https://github.com/sylviamic/python_bioinformagicks && \
	cd python_bioinformagicks && \
	make install && \
	make clean

ENTRYPOINT ["python"]