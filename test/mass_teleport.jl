# 서버에 있는 이미지 파일의 디렉터리
from = "//155.230.152.208/Immune_network/NI_article_fig"
# 로컬에 이미지 파일을 풀어 헤쳐서 저장할 디렉터리
to = "C:/test/immune/NI_article_fig"

# 서버 순회, 첫번째 원소는 제외해야 해서 [2:end]로 슬라이싱
walkeddir = collect(walkdir(from))[2:end]
for wd = walkeddir
    src = first(wd) .* "/" .* last(wd)
    # __을 /로 바꾸면 원래의 디렉터리 구조로 복원됨
    dst = replace.(src, from => to, "/fig" => "__fig")
    cp.(src, dst, force = true)
end
