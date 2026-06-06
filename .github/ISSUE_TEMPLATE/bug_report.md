---
name: Bug report
about: Report a problem so it can be reproduced and fixed
title: ''
labels: bug
assignees: joonho112
---

**Describe the bug**

A clear and concise description of what went wrong.

**Reproducible example**

Please include a minimal [reprex](https://reprex.tidyverse.org) — the smallest
self-contained code that reproduces the problem. Use the bundled example data or
the tiny fixture where possible so the report runs without a Gurobi license.

```r
library(gradepath)
# minimal code that reproduces the issue
```

**Expected behavior**

What you expected to happen instead.

**Solver / backend**

Which backend were you using (`gurobi` default, or `highs`)? Was Gurobi licensed
and on `PATH`?

**Session info**

<details>

```r
sessionInfo()
# or, if installed: sessioninfo::session_info()
```

</details>
