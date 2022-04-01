using PlotlyJS, LaTeXStrings


trace1 = scatter(

    x=[1, 2, 3, 4],

    y=[1, 4, 9, 16],

    mode="markers",

    name=L"$\alpha_{1c} = 352 \pm 11 \text{ km s}^{-1}$"

)


trace2 = scatter(

    x=[1, 2, 3, 4],

    y=[0.5, 2, 4.5, 8],

    name=L"$\beta_{1c} = 25 \pm 11 \text{ km s}^{-1}$"

)

layout = Layout(

    xaxis_title=L"$\sqrt{(n_\text{c}(t|{T_\text{early}}))}$",

    yaxis_title=L"$d, r \text{ (solar radius)}$"

)


plot([trace1, trace2], layout)

using Plots, PlotlyJS
plotlyjs()
PlotlyJS.plot(Plots.fakedata(50, 5), w = 3)

PlotlyJS.plot(1:4, [[1,4,9,16]*10000, [0.5, 2, 4.5, 8]],
           labels = [L"\alpha_{1c} = 352 \pm 11 \text{ km s}^{-1}";
                     L"\beta_{1c} = 25 \pm 11 \text{ km s}^{-1}"] |> permutedims,
           xlabel = L"\sqrt{(n_\text{c}(t|{T_\text{early}}))}",
           ylabel = L"d, r \text{ (solar radius)}",
           yformatter = :plain,
           extra_plot_kwargs = KW(
               :include_mathjax => "cdn",
               :yaxis => KW(:automargin => true),
               :xaxis => KW(:domain => "auto")
               ),
       )