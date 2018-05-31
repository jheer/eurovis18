const subvis = (function() {
  let width,
      height,
      data,
      svg,
      axisLayer,
      dataLayer,
      sim,
      radius = 10;

  const tableau10 = [
    "#4c78a8", "#f58518", "#e45756", "#72b7b2", "#54a24b",
    "#eeca3b", "#b279a2", "#ff9da6", "#9d755d", "#bab0ac"
  ];

  const color = d3.scaleOrdinal()
    .domain(["Accept", "Fast-Track", "Reject"])
    .range(tableau10);

  const typeAxis = d3.scaleBand()
    .domain(["Algorithm", "Design Study", "Evaluation", "System", "Model"]);

  const scoreAxis = d3.scaleLinear()
    .domain([1, 5]);

  function initialize(node, w, h) {
    width = w;
    height = h;

    svg = d3.select(node).append('svg')
      .attr('viewBox', `0 0 ${width} ${height}`)
      .attr('width', width)
      .attr('height', height);

    axisLayer = svg.append('g');
    dataLayer = svg.append('g');

    typeAxis.range([~~(0.15*height), ~~(0.85*height)]);
    scoreAxis.range([~~(0.1*width), ~~(0.9*width)]);

    d3.json('../data/submissions.json')
      .then(function(d) {
        data = d;
        window.parent.data = d;
        dataLayer.selectAll('circle')
          .data(data)
          .enter().append('circle')
            .style('fill', tableau10[0])
            .attr('cx', width / 2)
            .attr('cy', height / 2)
            .attr('r', radius);
        toBubble();
      })
      .catch(e => console.error(e));
  }

  function toBubble() {
    clearAxes();
    addTitle('162 Full Paper Submissions');

    dataLayer.selectAll('circle')
      .transition()
      .duration(1000)
      .style('fill', tableau10[0]);

    data.forEach(d => {
      d.x = width / 2;
      d.y = height / 2;
      d.vx = 0;
      d.vy = 0;
    });

    sim = d3.forceSimulation(data)
      .force('collide', d3.forceCollide(radius + 5))
      .force('x', d3.forceX(~~(width/2)))
      .force('y', d3.forceY(~~(height/2)));

    sim.on('tick', function() {
      dataLayer.selectAll('circle')
        .attr('cx', d => d.x)
        .attr('cy', d => d.y);
    });
  }

  function colorDecision() {
    clearAxes('g.title');
    addTitle('47 Accepted Papers (29%), 2 Fast Tracked to CGF');
    dataLayer.selectAll('circle')
      .transition()
      .duration(1000)
      .style('fill', d => color(d.decision));
  }

  function toScores() {
    clearAxes();
    addTitle('The Review Process...');

    let baseline = ~~(0.85 * height);

    axisLayer.append('g')
      .attr('transform', `translate(0, ${baseline})`)
      .style('opacity', 0)
      .call(d3.axisBottom(scoreAxis))
      .call(axisConfig)
    .transition()
      .duration(1000)
      .style('opacity', 1);

    baseline -= radius - 1;
    calculateStacks(data, 'score');

    dataLayer.selectAll('circle')
      .transition()
        .duration(2000)
        .delay((d, i) => 200 * d.score)
        .attr('cx', d => scoreAxis(d.score))
        .attr('cy', baseline)
        .style('fill', d => tableau10[0])
      .transition()
        .duration(1000)
        .attr('cy', d => baseline - (d.stack * 2 * radius));
  }

  function toTypes() {
    clearAxes();
    addTitle('Paper Types');

    let offset = ~~(0.1 * width);

    axisLayer.append('g')
      .attr('transform', `translate(${offset}, 0)`)
        .style('opacity', 0)
        .call(d3.axisLeft(typeAxis))
        .call(axisConfig)
      .transition()
        .duration(1000)
        .style('opacity', 1);

    offset += radius + 1;
    calculateStacks(data, 'type');

    dataLayer.selectAll('circle')
      .transition()
        .duration(1500)
        .delay((d, i) => i * 2)
        .attr('r', radius)
        .attr('cx', d => offset + (d.stack * 2 * radius))
        .attr('cy', d => typeAxis(d.type) + 0.5 * typeAxis.bandwidth());
  }

  function clearAxes(selector) {
    if (sim) {
      sim.stop();
      sim = null;
    }
    axisLayer.selectAll(selector || 'g')
      .transition()
      .style('opacity', 0)
      .remove();
  }

  function axisConfig(g) {
    g.selectAll('text')
      .style('font-size', '14px');
  }

  function addTitle(title) {
    const x = ~~(0.5 * width),
          y = ~~(0.1 * height) - 5;

    axisLayer.append('g')
      .attr('class', 'title')
      .attr('transform', `translate(${x}, ${y})`)
      .append('text')
      .text(title)
      .style('font-weight', 500)
      .style('font-size', '36px')
      .style('font-weight', 'bold')
      .style('text-anchor', 'middle')
      .style('opacity', 0)
    .transition()
      .duration(1000)
      .style('opacity', 1);
  }

  function calculateStacks(data, field) {
    let p = {};

    data.forEach(function(d) {
      const n = p[d[field]] || 0;
      p[d[field]] = (d.stack = n) + 1;
    });
  }

  return {
    initialize:    initialize,
    toBubble:      toBubble,
    toScores:      toScores,
    toTypes:       toTypes,
    toDecision:    colorDecision
  };

})();
