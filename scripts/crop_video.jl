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
