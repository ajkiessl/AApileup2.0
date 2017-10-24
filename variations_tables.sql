use testing;

drop table if exists substitutions;
drop table if exists sites;
drop table if exists templates;


create table templates (
    id serial primary key,
    protein text not null,
    length bigint unsigned not null
);

create table sites (
    id serial primary key,
    template_id bigint unsigned not null,
    position bigint unsigned not null,
    wild_type_AA char(1) not null,
    foreign key (template_id) references templates (id)
);

create table substitutions (
	id serial primary key,
	site_id bigint unsigned not null,
	substitution char(1) not null,
	count bigint unsigned not null,
    foreign key (site_id) references sites (id),
    unique key (site_id, substitution)
);

create or replace view AApileup as
    select templates.id as template_id,
       templates.protein as protein,
       sites.id as site_id,
       sites.position as position,
       sites.wild_type_AA as wild_type_AA,
       substitutions.substitution as substitution,
       substitutions.count as count
   from templates
       inner join sites on templates.id = sites.template_id
       inner join substitutions on sites.id = substitutions.site_id;