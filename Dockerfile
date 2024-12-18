FROM python:3.12.8

RUN \
	git clone https://github.com/sylviamic/python_bioinformagicks && \
	cd python_bioinformagicks && \
	make install && \
	make clean

ENTRYPOINT ["python"]
