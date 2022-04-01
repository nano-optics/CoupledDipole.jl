using WebIO, JSExpr

w = Scope(imports=["//cdnjs.cloudflare.com/ajax/libs/p5.js/0.5.11/p5.js"])
onmount(w, @js function (p5)
    function sketch(s)
        s.setup = () -> s.createCanvas(640, 480)

        s.draw = function ()
          if s.mouseIsPressed
            s.fill(0)
          else
            s.fill(255)
          end
          s.ellipse(s.mouseX, s.mouseY, 80, 80)
        end
    end
    @new p5(sketch, this.dom.querySelector(".container"))
end)
w(dom"div.container"())