using PackageCompiler

Pkg.activate(".")

import StepsProcessor

PackageCompiler.create_sysimage(:StepsProcessor;
                                sysimage_path="ExampleSysimagePrecompile.so",
                                precompile_execution_file="precompile_steps_processor.jl")
