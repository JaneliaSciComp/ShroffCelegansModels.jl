using VideoIO
using ColorTypes
using Makie: VideoStream
using Base64

# Makie-saved videos report a degenerate framerate (1//0) that libx264 rejects at
# codec-open ("Could not open codec: Return code -22"). Return a positive, finite
# framerate, falling back to `default` for any non-finite/non-positive value.
function sane_framerate(fr; default = 24, source = nothing)
    if isfinite(float(fr)) && float(fr) > 0
        return fr
    end
    @warn "Invalid video framerate $(fr); defaulting to $default" source
    return default
end

function crop_bounds(img; offset=50)
    d1 = sum(@view(img[(offset+1):end,(offset+1):end]), dims=2)
    d2 = sum(@view(img[(offset+1):end,(offset+1):end]), dims=1)
    ranges =
    (findfirst(!=(RGB(0,0,0)), d1)[1]:findlast(!=(RGB(0,0,0)), d1)[1]) .+ offset,
    (findfirst(!=(RGB(0,0,0)), d2)[2]:findlast(!=(RGB(0,0,0)), d2)[2]) .+ offset 
    ranges = map(ranges) do r
        if isodd(length(r))
            first(r):last(r)+1
        else
            r
        end
    end
end
function crop_video(
    filename::String,
    out_filename::String = replace(filename, ".mp4" => "_cropped.mp4");
    offset = 50,
    framerate = nothing,
    crf = 10,
)
    vio = openvideo(filename)
    N = counttotalframes(vio)
    skipframes(vio, N-1)
    last_frame = read(vio)
    bounds = crop_bounds(last_frame; offset)
    seekstart(vio)
    # Prefer an explicitly-passed framerate (e.g. vs.options.framerate); fall back
    # to the file's rate, sanitized (Makie files report a degenerate 1//0).
    _framerate = sane_framerate(isnothing(framerate) ? VideoIO.framerate(vio) : framerate; source = filename)
    # ffmpeg -i 2024_10_11_edited_xz_v2_cropped.mp4 -profile:v high422 -crf 17 -preset slow -c:v libx264 -pix_fmt yuv420p -an 2024_10_11_edited_xz_v5_cropped.mp4
    # Lower crf = sharper (less compression); 0 is lossless, 23 is the x264 default.
    # profile "high" matches the yuv420p we encode for broad browser playback
    # ("high422" requires 4:2:2 chroma; high422 + yuv420p fails codec-open, EINVAL -22).
    open_video_out(
        out_filename,
        @view(last_frame[bounds...]);
        codec_name = "libx264",
        encoder_options = (; crf, preset="slow", profile="high"),
        target_pix_fmt = VideoIO.AV_PIX_FMT_YUV420P,
        framerate = _framerate
    ) do writer
        for i in 1:N
            write(writer, @view(read(vio)[bounds...]))
        end
    end
    close(vio)
end

function crop_video(vs::VideoStream)
    # From https://github.com/MakieOrg/Makie.jl/blob/2acd423116b3c28b0d51762be2747d7ebb3eee84/src/recording.jl#L182C1-L193C1
    # The MIT License (MIT) Copyright (c) 2018-2021: Simon Danisch, Julius Krumbiegel.
    mktempdir() do dir
        raw_path = save(joinpath(dir, "video.mp4"), vs)
        path = joinpath(dir, "cropped_video.mp4")
        # the saved file's framerate metadata is degenerate (1//0); pass the
        # VideoStream's real framerate so the re-encode uses a valid rate.
        crop_video(raw_path, path; framerate = vs.options.framerate)

        # <video> only supports infinite looping, so we loop forever even when a finite number is requested
        loopoption = vs.options.loop ≥ 0 ? (;loop=true) : (;)
        source = DOM.source(; src="data:video/x-m4v;base64,$(base64encode(open(read,path)))", type="video/mp4")
        # TODO: use loop option
        return DOM.video(source; autoplay=true, controls=true, loopoption...);
    end
end
