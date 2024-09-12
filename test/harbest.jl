using CSV, DataFrames, Harbest, Snowball

❎ = ["", "Games & Quizzes", "Shop", "Books", "", "Merch", "Settings", "My Words", "Recents", "Account", "More", "Join MWU", "Username"]

ab = ['a':'z'; '0']
pn = [53,23,55,30,30,18,17,35,27,2,6,23,44,19,19,64,2,21,50,28,7,13,5,2,1,3,1]

term = String[]
for (i, x) = enumerate(ab)
    for y = 1:pn[i]
        page = read_html("https://www.merriam-webster.com/browse/medical/$x/$y")
        word = setdiff(html_elements(page, ["span"]) |> html_text3, ❎)
        append!(term, word)
        print("$x/$y-")
    end
end
term = [stem_all(Stemmer("english"), trm) for trm in lowercase.(term)]
unique!(term)
gowords = DataFrame(; term)
CSV.write("gowords.csv", gowords, bom = true)
