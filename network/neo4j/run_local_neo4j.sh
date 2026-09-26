#!/usr/bin/env bash
# =====================================================================================================
# network/neo4j/run_local_neo4j.sh — START A LOCAL NEO4J (DOCKER) AND LOAD THE EXERCISE NETWORKS
# =====================================================================================================
# PURPOSE: one command from pipeline outputs to a browsable graph at http://localhost:7474, for building and
#   testing the visualiser. Uses the official neo4j:5-community image; no password (local use only).
# HOW TO RUN (from the repo root, after the pipeline):   bash network/neo4j/run_local_neo4j.sh
#   Stop and remove:   docker rm -f hackathon-neo4j        (the database lives in the container only)
# INPUTS: the CSVs written by export_neo4j.R ($NEO4J_IMPORT, default $HACK_OUT/neo4j_import)
# =====================================================================================================
set -euo pipefail

# Folder of this script (so it works from anywhere) and the export folder.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HACK_OUT="${HACK_OUT:-$HOME/Desktop/output/hackathon-2026-track1/network}"
NEO4J_IMPORT="${NEO4J_IMPORT:-$HACK_OUT/neo4j_import}"
NAME="hackathon-neo4j"

# 1. Refresh the export from the current pipeline outputs.
Rscript "$HERE/export_neo4j.R"

# 2. Start a fresh container with the export mounted as Neo4j's import folder (read-only).
docker rm -f "$NAME" >/dev/null 2>&1 || true
docker run -d --name "$NAME" -p 7474:7474 -p 7687:7687 -e NEO4J_AUTH=none \
  -v "$NEO4J_IMPORT":/import:ro -v "$HERE":/scripts:ro neo4j:5-community >/dev/null

# 3. Wait until the database answers (up to 2 minutes), stopping early if the container dies.
for i in $(seq 1 60); do
  if docker exec "$NAME" cypher-shell "RETURN 1" >/dev/null 2>&1; then break; fi
  docker ps -q -f name="$NAME" | grep -q . || { echo "Neo4j container exited:"; docker logs "$NAME" | tail -20; exit 1; }
  if [ "$i" -eq 60 ]; then echo "Neo4j did not start within 2 minutes"; exit 1; fi
  sleep 2
done

# 4. Load the graph and print the counts (compare with README "Expected counts").
docker exec "$NAME" cypher-shell -f /scripts/import.cypher
echo "Neo4j Browser: http://localhost:7474  (no login; try the queries in network/neo4j/queries.cypher)"
