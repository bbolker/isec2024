talk: shapeconst_talk.rmd waterbug.Rout reedfrog.Rout
	Rscript -e "rmarkdown::render('shapeconst_talk.rmd')"

%.Rout: %.R
	R CMD BATCH --vanilla $<

%.html: %.rmd
	Rscript -e "rmarkdown::render('$<')"

clean:
	rm *~ \#*
