# Config Change Summary

## 2024-09-05 v1 → 2026-03-19 v2

Settings (folderpaths, voxel scale, outlier removal, interpolation, warping, smoothing, mipav_output) are unchanged.

### Removed paths (5)

**JCC596_NU** — 3 old paths replaced by re-processed `Untwisting_Redo` equivalents:
- `JCC596_NU/Untwisting/082619_JCC596_NU/JCC596_NU/Pos3/SPIMB_result/Reg_Sample/Decon_registered`
- `JCC596_NU/Untwisting/091119_JCC596_NU/JCC596_NU/JCC596_NU/JCC596_NU/Pos2/SPIMB_result/Reg_Sample/Decon_registered`
- `JCC596_NU/Untwisting/091119_JCC596_NU/JCC596_NU/JCC596_NU/JCC596_NU/Pos3/SPIMB_result/Reg_Sample/Decon_registered`

**DCR4221** — UTP 1 dataset dropped; strain now has 4 datasets (UTP 2–5):
- `Non-model C elegans folders/Untwisting_Paper/DCR4221/063014_lattices_Javier - UTP 1/Decon_reg`

**RW10752_NU Pos2** — old path replaced:
- `RW10752_NU/Untwisting/031219_RW10752_NU/RW10752_NU/RW10752_NU/Pos2/Old/Decon_registered`

### Added paths (13)

**JCC596_NU** — 3 replacement paths under `Untwisting_Redo`:
- `JCC596_NU/Untwisting_Redo/082619_Pos3`
- `JCC596_NU/Untwisting_Redo/091119_Pos2`
- `JCC596_NU/Untwisting_Redo/091119_Pos3`

**RW10752_NU Pos2** — replacement path:
- `RW10752_NU/Untwisting/031219_RW10752_NU/RW10752_NU/RW10752_NU/Pos2/SPIMB/For_Tracking/For_Tracking`

**RW10131_retracked** — new strain, 3 positions (2024 SLS268 imaging):
- `RW10131/Data/2024_SLS268/20240401/RW10131_SLS6_New/Pos4/SPIMB/Reg_Sample/For_Tracking`
- `RW10131/Data/2024_SLS268/20240429/SLS268_RW10131_SLS6/Pos1/SPIMB/Reg_Sample/For_Tracking`
- `RW10131/Data/2024_SLS268/20240507/Pos1/SPIMA/Reg_Sample/For_Tracking`

**RW10598_retracked** — new strain, 3 positions (2023 re-tracking):
- `RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos1/SPIMB/Reg_Sample/For_Tracking`
- `RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos2/SPIMB/Reg_Sample/For_Tracking`
- `RW10598/2023_Data/Tracking/20230719/RW10598_NU/Pos4/SPIMB/Reg_Sample/For_Tracking`

**RW10896_retracked** — new strain, 3 positions (2023 imaging):
- `RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos1/SPIMB/Registered_Volumes/For_Tracking`
- `RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos2/SPIMB/Registered_Volumes/For_Tracking`
- `RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos3/SPIMB/Registered_Volumes/For_Tracking`
