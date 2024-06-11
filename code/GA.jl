using Plots

obstaculo(x, w, y, h) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

function dentro(x, y, obstaculo)
	return x ≥ obstaculo.x[1] && x ≤ obstaculo.x[2] && y ≥ obstaculo.y[1] && y ≤ obstaculo.y[3]
end

function criar_ecra(obstaculos, largura; Δ = 1)
    E = Matrix{Bool}(undef, (largura / Δ + 1) |> floor |> Int, (largura / Δ + 1) |> floor |> Int)
    for i in eachindex(IndexCartesian(), E)
        E[i] = dentro.((i[2]) * Δ, (i[1]) * Δ, obstaculos) |> !any
    end
    return E
end

function coordenadas_ecra(E, Δ, tipo)
    x, y = [], []
    for i in eachindex(IndexCartesian(), E)
        if E[i] == tipo
            push!(x, (i[2] - 1) * Δ + 1)
            push!(y, (i[1] - 1) * Δ + 1)
        end
    end
    return x, y
end

function mapear(obstaculos, largura; ecra=[], Δ = 1)
	p = plot(xlim=[1, largura + 1], ylim=[1, largura + 1], leg=:false)
	p = scatter!([largura], [largura], alpha=0.0)
	for obstaculo ∈ obstaculos
		p = plot!(obstaculo, color=:black)
	end
	if !isempty(ecra)
		p = scatter!(coordenadas_ecra(ecra, Δ, false)..., color=:red)
		p = scatter!(coordenadas_ecra(ecra, Δ, true)..., color=:blue)
	end
	return p
end

@enum Movimentos MOVIMENTO_TIPO_LINHA=0 MOVIMENTO_TIPO_COLUNA=1
@enum Ordens HORIZONTAL_VERTICAL=0 VERTICAL_HORIZONTAL=1

criar_cromossomo(n) = [rand(instances(Movimentos)); [1]; rand(1:n, n-2); [n]; rand(instances(Ordens), n - 1); rand(1:n, 2)]

function seccionar_cromossomo(cromossomo)
    n = Int((length(cromossomo) / 2) - 1)
    TIPO = 1
    ROTA = 2:n+1
    ORDEM = n+2:2*n
    TROCA = 2*n+1:2*n+2
    return cromossomo[TIPO], cromossomo[ROTA], cromossomo[ORDEM], cromossomo[TROCA]
end

function criar_populacao(tamanho, npontos)
    populacao = []
    for i in 1:tamanho
        push!(populacao, criar_cromossomo(npontos))
    end
    return populacao
end

function mapear_cromossomo(grafico, cromossomo, Δ, cor)
    tipo, rota, ordem, troca = seccionar_cromossomo(cromossomo)
    anterior = (0, 0)
    for (i, ponto) ∈ enumerate(rota)
        ax, ay = anterior
        if i ∈ troca
            if tipo == MOVIMENTO_TIPO_LINHA
                tipo = MOVIMENTO_TIPO_COLUNA
            else
                tipo = MOVIMENTO_TIPO_LINHA
            end
        end
        novo = tipo == MOVIMENTO_TIPO_LINHA ? (ponto, i) : (i, ponto)
        nx, ny = novo
        if anterior ≠ (0, 0)
            if ordem[i - 1] == HORIZONTAL_VERTICAL
                grafico = plot(grafico, Δ * ([ax; nx]), Δ * ([ay; ay]), color=cor, linewidth=5)
                grafico = plot(grafico, Δ * ([nx; nx]), Δ * ([ay; ny]), color=cor,  linewidth=5)
            else
                grafico = plot(grafico, Δ * ([ax; ax]), Δ * ([ay; ny]), color=cor,  linewidth=5)
                grafico = plot(grafico, Δ * ([ax; nx]), Δ * ([ny; ny]), color=cor, linewidth=5)
            end
        end
        anterior = novo
    end
    return grafico
end

function trocar_tipo_movimento(tipo)
    if tipo == MOVIMENTO_TIPO_LINHA
        tipo = MOVIMENTO_TIPO_COLUNA
    else
        tipo = MOVIMENTO_TIPO_LINHA
    end
    return tipo
end

function gerar_coordenadas_iniciais(tipo_movimento, rota)
    local x, y
    if tipo_movimento == MOVIMENTO_TIPO_LINHA
        x, y = first(rota), 1
    else
        x, y = 1, first(rota)
    end
    return x, y
end

function gerar_proximas_coordenadas(tipo_movimento, indice_gene, ponto)
    local nx, ny
    if tipo_movimento == MOVIMENTO_TIPO_LINHA
        nx, ny = ponto, indice_gene + 1
    else
        nx, ny = indice_gene + 1, ponto
    end
    return nx, ny
end

function caminho(x, nx)
    if sign(nx - x) == 0
        return x:nx
    end
    return x:sign(nx-x):nx
end

function fitness_parametros(cromossomos, ecra)
    sentido, direcao, distancias, colisoes, curvas = nothing, nothing, [], [], []
    for cromossomo ∈ cromossomos
        tipo, rota, ordem, troca = seccionar_cromossomo(cromossomo)
        x, y = gerar_coordenadas_iniciais(tipo, rota)
        dist, col, cur = 0, ecra[x,y] ? 0 : 1, 0
        for (i, ponto) ∈ enumerate(rota[2:end])
            if (i + 1) ∈ troca
                tipo = trocar_tipo_movimento(tipo)
            end
            nx, ny = gerar_proximas_coordenadas(tipo, i, ponto)
            dist += abs(nx - x + ny - y)
            if ordem[i] == HORIZONTAL_VERTICAL
                direcao = :horizontal
            else
                direcao = :vertical
            end
            while x ≠ nx || y ≠ ny
                if direcao == :horizontal
                    sinal = sign(nx - x)
                    if sinal == 0
                        direcao = :vertical
                        continue
                    end
                    xant, x = x, x + sinal
                    col += !ecra[y, x] && ecra[y, xant] ? 1 : 0
                    if sinal == 1 && sentido ≠ :leste
                        sentido = :leste
                        cur += isnothing(sentido) ? 0 : 1
                    elseif sinal == -1 && sentido ≠ :oeste
                        sentido = :oeste
                        cur += isnothing(sentido) ? 0 : 1
                    end
                end
                if direcao == :vertical
                    sinal = sign(ny - y)
                    if sinal == 0
                        direcao = :horizontal
                        continue
                    end
                    yant, y = y, y + sinal
                    col += !ecra[y, x] && ecra[yant, x] ? 1 : 0
                    if sinal == 1 && sentido ≠ :norte
                        cur += isnothing(sentido) ? 0 : 1
                        sentido = :norte
                    elseif sinal == -1 && sentido ≠ :sul
                        cur += isnothing(sentido) ? 0 : 1
                        sentido = :sul
                    end
                end
            end
        end
        push!(distancias, dist)
        push!(colisoes, col)
        push!(curvas, cur)
    end
    return distancias, colisoes, curvas
end

function mapear_intervalo!(valores, imin, imax, vmin, vmax)
    coef = (imax - imin) / (vmax - vmin)
    for i ∈ eachindex(valores)
        valores[i] = coef * (valores[i] - vmin) + imin
    end
end

function normalizar_parametros!(distancias, colisoes, curvas)
    mapear_intervalo!(distancias, 1, 0, minimum(distancias), maximum(distancias))
    mapear_intervalo!(colisoes, 1, 0, minimum(colisoes), maximum(colisoes))
    mapear_intervalo!(curvas, 1, 0, minimum(curvas), maximum(curvas))
end

function fitness(populacao, ecra; L=1, T=2)
	F = []
	distancias, colisoes, curvas = fitness_parametros(populacao, ecra)
    col = copy(colisoes)
	normalizar_parametros!(distancias, colisoes, curvas)
	for i in eachindex(populacao)
        fit = colisoes[i] * (L * distancias[i] + T * curvas[i]) * 100 / (L + T)
        if col[i] > 0
            fit *= 0.1 / (col[i]^2)
        end
		push!(F, fit)
	end
    ordem = sortperm(F)
    F = F[ordem]
    populacao[1:end] = populacao[ordem]
	return F
end

function criar_intervalo(populacao, rankings)
    soma_rankings = sum(rankings)
    soma_cumulativa = 0
    intervalos = [0.]
    for ranking ∈ rankings
        push!(intervalos, soma_cumulativa + ranking / soma_rankings)
        soma_cumulativa = last(intervalos)
    end
    rankings[end] = 1.0
    return intervalos
end

function escolher_pai(intervalos)
    selecao = rand()
    for i ∈ 2:length(intervalos)
        if intervalos[i - 1] ≤ selecao ≤ intervalos[i]
            selecao = i - 1
            break
        end
    end
    return selecao
end

function selecionar_casais(populacao)
    casais = []
    n = length(populacao)
    soma_rankings = n*(n+1)/2
    intervalos = criar_intervalo(populacao, collect(1:n))
    for i ∈ 1:length(populacao)/2
        pai1 = escolher_pai(intervalos)
        rankings = collect(1:n)
        rankings[pai1] = 0
        intervalos2 = criar_intervalo(populacao, rankings)
        pai2 = escolher_pai(intervalos2)
        push!(casais, (populacao[pai1], populacao[pai2]))
    end
    return casais
end

function crossover(casais)
    filhos = []
    for casal ∈ casais
        ponto_crossover = rand(3:(length(casal[1])-4))
        push!(filhos, [casal[1][1:ponto_crossover];casal[2][ponto_crossover+1:end]])
        push!(filhos, [casal[2][1:ponto_crossover];casal[1][ponto_crossover+1:end]])
    end
    return filhos
end

function mutacao(populacao, probabilidade)
    n = length(populacao[1])
    N = Int(n / 2 - 1)
    mutacoes = 0
    for cromossomo ∈ populacao
        if rand() ≤ probabilidade
            mutacoes += 1
            if cromossomo[1] == MOVIMENTO_TIPO_COLUNA
                cromossomo[1] = MOVIMENTO_TIPO_LINHA
            else
                cromossomo[1] = MOVIMENTO_TIPO_COLUNA
            end
        end
        for i ∈ 3:N-1
            if rand() ≤ probabilidade
                mutacoes += 1
                cromossomo[i] = rand(1:N)
            end
        end
        for i ∈ N+2:n-2
            if rand() ≤ probabilidade
                mutacoes += 1
                if cromossomo[i] == HORIZONTAL_VERTICAL
                    cromossomo[i] = VERTICAL_HORIZONTAL
                else
                    cromossomo[i] = HORIZONTAL_VERTICAL
                end
            end
        end
        for i ∈ n-1:n
            if rand() ≤ probabilidade
                mutacoes += 1
                cromossomo[i] = rand(1:N)
            end
        end
    end
    println("Mutações aplicadas a $(100 * mutacoes / (length(populacao[1]) * length(populacao)))% da geração")
end

function algoritmo_genetico(tamanho_populacao, ecra, probabilidade_mutacao, iteracoes, intervalo, grafico, Δ, cor)
    populacao = criar_populacao(tamanho_populacao, size(ecra, 1))
    F = fitness(populacao, ecra)
    animacao = Animation()
    for i ∈ 1:iteracoes
        println("Média fitness da geração $i = $(sum(F) / length(populacao))")
        casais = selecionar_casais(populacao)
        populacao = crossover(casais)
        mutacao(populacao, probabilidade_mutacao)
        F = fitness(populacao, ecra)
        if i % intervalo == 0
            p = mapear_cromossomo(grafico, populacao[end], Δ, cor)
            title!(p, "Geração $i")
            frame(animacao, p)
        end
    end
    gif(animacao, "animacao.gif", fps=1)
    return populacao
end

##

# obs = [obstaculo(0, 5, 6, 4),
# 		obstaculo(4, 3, 0, 4),
# 		obstaculo(8, 1, 3, 7)]
# largura = 10
# Δ = 1

# obs1 = [obstaculo(0, 5, 6, 4),
# 		obstaculo(4, 3, 0, 4),
# 		obstaculo(8, 1, 3, 7)]
# origem1 = [1,1]
# destino1 = [10,10]
# largura1 = 10
# altura1 = 10
# Δexemplo = 0.75

# obs1 = [obstaculo(1, 1, 0, 4),
# 		obstaculo(3, 1, 2, 5),
# 		obstaculo(3, 3, 2, 1),
#          obstaculo(5, 2, 4, 1)]
# origem1 = [1,1]
# destino1 = [7,7]
# largura1 = 7
# altura1 = 7
# Δexemplo = 0.25

obs = [obstaculo(1.5, 2, 1, 3)]
largura = 4
∆ = 1
E = criar_ecra(obs, largura, Δ=Δ)
grafico = mapear(obs, largura, ecra=E, Δ=Δ)
graf = mapear_cromossomo(grafico, [MOVIMENTO_TIPO_COLUNA, 1, 5, 5, 5, 5, 1, 1, HORIZONTAL_VERTICAL, HORIZONTAL_VERTICAL, 0, 0], Δ, :green)
graf = plot(graf, [1, 2, 5], [1, 5, 5], color=:maroon1, linewidth = 5)
# E = criar_ecra(obs, largura, Δ = Δ)
# grafico = mapear(obs, largura, ecra=E, Δ = Δ)
# graf = mapear_cromossomo(grafico, [MOVIMENTO_TIPO_LINHA, 1, 3, 1, 2, 2, 6, 0, 0, 0, 0, 0, 0, 0], 1, :green)
# populacao = algoritmo_genetico(50, Eexemplo, 0.05, 10000, 100, grafico, Δexemplo, :green)
# mapear_cromossomo(grafico, populacao[20], Δexemplo, :green)