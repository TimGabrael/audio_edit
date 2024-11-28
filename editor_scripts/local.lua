
local function close_windows(name)
  -- Get a list of all windows
  local cmd = 'taskkill /IM ' .. name .. ' /F'
  local results = vim.fn.system(cmd)
  print(results)
end

vim.keymap.set('n', '<F7>', 
function()
    close_windows("audio_edit.exe")
    local root_dir = vim.fn.getcwd()
    local game_cmd1 = 'start ' .. root_dir .. '/bin/audio_edit.exe '
    vim.fn.system(game_cmd1)
end, { noremap = true, silent = true})

vim.keymap.set('n', '<F5>', 
function()
    close_windows("audio_edit.exe")
    local root_dir = vim.fn.getcwd()
    local gpt_cmd = 'start cmd /c "(cmake --build ' .. root_dir .. '/build/ > build_output.txt 2>&1) & findstr /i /c:\"error\" build_output.txt && (pause) || del build_output.txt"'
    local compile_cmd = 'start cmake --build ' .. root_dir .. '/build/'
    vim.fn.system(gpt_cmd)
end, { noremap = true, silent = true})

