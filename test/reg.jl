corpus = read("BTKAeLqLMw.txt", String)

reg = r"\d+ [A-Z]{7,}"

section = collect(eachmatch(reg, corpus))[1:end-1]
position = findall(reg, corpus)
volume = first.(position)[2:end] - last.(position)[1:end-1]
name = getproperty.(section, :match)
["'$n' 섹션의 글자 수: $v" for (n, v) in zip(name, volume)]

