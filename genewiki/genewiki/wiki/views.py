from django.template import RequestContext
from django.shortcuts import get_object_or_404, render_to_response, redirect
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.views.decorators.http import require_http_methods
from django.http import HttpResponse

from genewiki.mapping.models import Relationship

from genewiki.wiki.models import Article

from genewiki.wiki.textutils import create, interwiki_link

from datetime import datetime, timedelta


@require_http_methods(['GET', 'POST'])
def article_create(request, entrez_id):
    results = create(entrez_id)

    # We failed to gather information then return the ID error
    if results is None:
        return HttpResponse('Invalid or missing Entrez Identifier')

    titles = results.get('titles')
    article = Article.objects.get_infobox_for_entrez(entrez_id)

    vals = {'results': results,
            'article': article,
            'titles': titles,
            'title': Relationship.objects.get_title_for_entrez(entrez_id),
            'entrez': entrez_id,}

    if request.method == 'POST':
        # Only assign this 'title' var internally if the online article status is False (not a Wikipedia page)
        uploadopt = request.POST.get('page_type')
        if uploadopt is None:
            return HttpResponse('Must select title option.')

        title = titles[uploadopt][0] if titles[uploadopt][1] is False else None

        # The page title that they wanted to create is already online
        if title is None:
            return HttpResponse('Article or template already exists.')

        vals['title'] = title
        content = results['stub']
        Article.objects.get_or_create(title=title, text=content, article_type=Article.PAGE, force_update=True)

        # create corresponding talk page with appropriate project banners
        talk_title = 'Talk:{0}'.format(title)
        talk_content = """{{WikiProjectBannerShell|
                          {{WikiProject Gene Wiki|class=stub|importance=low}}
                          {{Wikiproject MCB|class=stub|importance=low}}
                          }}"""
        Article.objects.get_or_create(title=talk_title, text=talk_content, article_type=Article.TALK, force_update=True)
        #create interwiki link
        link = interwiki_link(entrez_id, title)
     

        # Save the entrez_id to title mapping for future reference
        if not Relationship.objects.filter(entrez_id=entrez_id).exists():
           Relationship.objects.create(entrez_id=entrez_id, title=title)

        return redirect('genewiki.wiki.views.article_create', entrez_id)

    return render_to_response('wiki/create.jade', vals, context_instance=RequestContext(request))

