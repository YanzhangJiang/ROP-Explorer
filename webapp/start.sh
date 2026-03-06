#!/bin/bash
# ROP Explorer Startup Script

cd "$(dirname "$0")"

PORT=${ROP_PORT:-8088}

echo "Starting ROP Explorer Web Server..."
echo "Server will be available at: http://localhost:$PORT"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

julia -t auto --project=. server.jl
