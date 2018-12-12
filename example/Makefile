

all: logexpr logexpr_steps zscores zscores_steps

logexpr:
	./__test_scripts/run.logexpr.R
	touch logexpr

logexpr_steps:
	./__test_scripts/run.logexpr.R.steps
	touch logexpr_steps

zscores:
	./__test_scripts/run.zscores.R
	touch zscores

zscores_steps:
	./__test_scripts/run.zscores.R.steps
	touch zscores_steps




clean:
	rm -f logexpr logexpr_steps zscores zscores_steps
	rm -rf ./__test_outdir.*
