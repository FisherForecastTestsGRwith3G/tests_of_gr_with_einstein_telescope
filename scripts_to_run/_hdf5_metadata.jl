"""
    write_top_level_output_metadata!(file; repo_root=joinpath(@__DIR__, ".."))

Attach repository-level provenance metadata to the root of a newly created
HDF5 file.

The stored attributes are:

- `date`: local calendar date when the file was created, formatted as
  `yyyy-mm-dd`
- `githash`: full commit hash of the repository HEAD at creation time

If the Git hash cannot be determined, `githash` is set to `"unknown"`.
"""
function write_top_level_output_metadata!(file; repo_root::AbstractString=joinpath(@__DIR__, ".."))
    attrs(file)["date"] = Dates.format(Dates.now(), dateformat"yyyy-mm-dd")

    githash = try
        readchomp(Cmd(`git rev-parse HEAD`, dir=repo_root))
    catch
        "unknown"
    end

    attrs(file)["githash"] = githash
    return nothing
end
