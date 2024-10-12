using VideoIO
using ColorTypes
using Makie: VideoStream
using Base64

function crop_bounds(img)
    d1 = sum(@view(img[51:end,51:end]), dims=2)
    d2 = sum(@view(img[51:end,51:end]), dims=1)
    ranges =
    (findfirst(!=(RGB(0,0,0)), d1)[1]:findlast(!=(RGB(0,0,0)), d1)[1]) .+ 50,
    (findfirst(!=(RGB(0,0,0)), d2)[2]:findlast(!=(RGB(0,0,0)), d2)[2]) .+ 50
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
    out_filename::String = replace(filename, ".mp4" => "_cropped.mp4")
)
    vio = openvideo(filename)
    N = counttotalframes(vio)
    skipframes(vio, N-1)
    last_frame = read(vio)
    bounds = crop_bounds(last_frame)
    seekstart(vio)
    # ffmpeg -i 2024_10_11_edited_xz_v2_cropped.mp4 -profile:v high422 -crf 17 -preset slow -c:v libx264 -pix_fmt yuv420p -an 2024_10_11_edited_xz_v5_cropped.mp4
    open_video_out(
        out_filename,
        @view(last_frame[bounds...]);
        codec_name = "libx264",
        encoder_options = (; crf=17, preset="slow", profile="high422"),
        target_pix_fmt = VideoIO.AV_PIX_FMT_YUV420P
    ) do writer
        for i in 1:N
            write(writer, @view(read(vio)[bounds...]))
        end
    end
end

function crop_video(vs::VideoStream)
    # From https://github.com/MakieOrg/Makie.jl/blob/2acd423116b3c28b0d51762be2747d7ebb3eee84/src/recording.jl#L182C1-L193C1
    # The MIT License (MIT) Copyright (c) 2018-2021: Simon Danisch, Julius Krumbiegel.
    mktempdir() do dir
        raw_path = save(joinpath(dir, "video.mp4"), vs)
        path = joinpath(dir, "cropped_video.mp4")
        crop_video(raw_path, path)

        # <video> only supports infinite looping, so we loop forever even when a finite number is requested
        loopoption = vs.options.loop ≥ 0 ? (;loop=true) : (;)
        source = DOM.source(; src="data:video/x-m4v;base64,$(base64encode(open(read,path)))", type="video/mp4")
        # TODO: use loop option
        return DOM.video(source; autoplay=true, controls=true, loopoption...);
    end
end
