all : test

clean:
	test -d index && rm -rf index
	test -d query && rm -rf query

test:
	pytest --git-aware --symlink
