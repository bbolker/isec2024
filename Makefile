talk:
	Rscript -e "rmarkdown::render('shapeconst_talk.rmd')"

%.Rout: %.R
	R CMD BATCH --vanilla $<