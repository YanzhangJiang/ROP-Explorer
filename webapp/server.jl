#!/usr/bin/env julia

import Pkg

# Make the server resilient to how it is launched.
# This allows `julia webapp/server.jl` as well as `--project=webapp`.
Pkg.activate(@__DIR__; io=devnull)

using ROPExplorerBackend

ROPExplorerBackend.main()
