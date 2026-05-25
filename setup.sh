#!/usr/bin/env bash
# Setup script — run once from inside /home/marko-b2/upwork_portfolio
# Usage: bash setup.sh YOUR_GITHUB_USERNAME

set -e

GH_USER="${1:-USERNAME}"
REPO_NAME="scrnaseq-portfolio"

echo "==> Setting up scRNA-seq portfolio for GitHub user: $GH_USER"

# 1. Git init
if [ ! -d .git ]; then
    git init -b main
    echo "    git initialized"
else
    echo "    git already initialized"
fi

# 2. Configure remote
if git remote get-url origin >/dev/null 2>&1; then
    echo "    remote 'origin' already set: $(git remote get-url origin)"
else
    git remote add origin "git@github.com:${GH_USER}/${REPO_NAME}.git"
    echo "    remote 'origin' set to git@github.com:${GH_USER}/${REPO_NAME}.git"
fi

# 3. First commit
git add .
if git diff --cached --quiet; then
    echo "    nothing to commit"
else
    git commit -m "Initial portfolio scaffold: 10 scRNA-seq projects, environment, templates"
    echo "    initial commit created"
fi

# 4. Conda env
echo ""
echo "==> Next steps:"
echo "    1. Create the GitHub repo (empty, no README):"
echo "       https://github.com/new  →  name: $REPO_NAME"
echo ""
echo "    2. Push:"
echo "       git push -u origin main"
echo ""
echo "    3. Create conda environment:"
echo "       mamba env create -f environment.yml"
echo "       mamba activate scportfolio"
echo ""
echo "    4. Register Jupyter kernels:"
echo "       python -m ipykernel install --user --name scportfolio --display-name 'Python (scportfolio)'"
echo "       Rscript -e 'IRkernel::installspec(name=\"r-scportfolio\", displayname=\"R (scportfolio)\")'"
echo ""
echo "    5. Start with Project 01 (PBMC 3K) as sanity check"
