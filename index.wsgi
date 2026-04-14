# -*- coding: utf-8 -*-
"""
ClusterChirp — Clean Production WSGI

Lightweight version for production use - no debugging overhead
"""

import os
import sys
import site

# ----------------------------
# CONFIGURATION
# ----------------------------
VENV_SITE_PACKAGES = "/hpc/users/rawalo01/heatmapBackendVirtual/lib/python3.9/site-packages"
PROJECT_DIR         = "/hpc/users/rawalo01/www/clusterchirp/backend"
DJANGO_SETTINGS     = "HeatmapBackend.settings"

# Runtime directory for matplotlib (avoid NFS issues)
RUNTIME_DIR = os.path.join(PROJECT_DIR, "_runtime")
MPL_TMP     = os.path.join(RUNTIME_DIR, "mpltmp")

# ----------------------------
# SETUP
# ----------------------------

# Configure matplotlib temp directory (lightweight)
try:
    os.makedirs(MPL_TMP, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", MPL_TMP)
except Exception:
    pass  # Non-critical, continue if it fails

# Add project and virtual environment to Python path
if PROJECT_DIR not in sys.path:
    sys.path.insert(0, PROJECT_DIR)
site.addsitedir(VENV_SITE_PACKAGES)

# Set Django settings module
os.environ.setdefault("DJANGO_SETTINGS_MODULE", DJANGO_SETTINGS)

# Load OpenAI API key from project-local file (gitignored, never committed)
try:
    key_path = os.path.join(PROJECT_DIR, ".openai_key")
    with open(key_path) as f:
        os.environ["OPENAI_API_KEY"] = f.read().strip()
except FileNotFoundError:
    pass  # AI chat feature will return a 503 if key is missing

# ----------------------------
# WSGI APPLICATION
# ----------------------------
from django.core.wsgi import get_wsgi_application
application = get_wsgi_application()
